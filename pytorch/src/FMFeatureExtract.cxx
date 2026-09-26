// FMFeatureExtract: frame -> FM input canvas -> ITensorForward -> per-pixel feature sidecar.
// See the header and wcfm/docs/01 sec 4.2, wcfm/docs/04.
#include "WireCellPytorch/FMFeatureExtract.h"
#include "WireCellPytorch/Torch.h"  // c10::Half

#include "WireCellAux/FrameTools.h"
#include "WireCellAux/PlaneTools.h"
#include "WireCellAux/SimpleTensor.h"
#include "WireCellAux/SimpleTensorSet.h"

#include "WireCellUtil/Array.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/NamedFactory.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <set>

WIRECELL_FACTORY(FMFeatureExtract, WireCell::Pytorch::FMFeatureExtract,
                 WireCell::INamed, WireCell::IFrameTensorSet, WireCell::IConfigurable)

using namespace WireCell;

Pytorch::FMFeatureExtract::FMFeatureExtract()
    : Aux::Logger("FMFeatureExtract", "torch")
{
}

Pytorch::FMFeatureExtract::~FMFeatureExtract() {}

WireCell::Configuration Pytorch::FMFeatureExtract::default_configuration() const
{
    Configuration cfg;
    cfg["anode"] = m_anode_tn;
    cfg["forward"] = m_forward_tn;
    cfg["input_tag"] = m_input_tag;
    cfg["planes"] = Json::arrayValue;
    for (int p : m_planes) cfg["planes"].append(p);
    cfg["tick0"] = m_tick0;
    cfg["nticks"] = m_nticks;
    cfg["tick_span"] = m_tick_span;
    cfg["input_scale"] = m_input_scale;
    cfg["input_offset"] = m_input_offset;
    cfg["active_threshold"] = m_active_threshold;
    cfg["view_norm"] = Json::arrayValue;
    for (const auto& vn : m_view_norm) {
        Configuration one = Json::arrayValue;
        one.append(vn[0]);
        one.append(vn[1]);
        cfg["view_norm"].append(one);
    }
    cfg["min_canvas"] = m_min_canvas;
    cfg["bbox_pad"] = m_bbox_pad;
    cfg["feature_dim"] = m_feature_dim;
    cfg["max_dense_pixels"] = (Json::Int64) m_max_dense_pixels;
    cfg["halo"] = m_halo;
    cfg["store_half"] = m_store_half;
    cfg["provenance"] = m_provenance;
    return cfg;
}

void Pytorch::FMFeatureExtract::configure(const WireCell::Configuration& cfg)
{
    m_anode_tn = get(cfg, "anode", m_anode_tn);
    m_forward_tn = get(cfg, "forward", m_forward_tn);
    m_input_tag = get(cfg, "input_tag", m_input_tag);
    if (cfg.isMember("planes")) {
        m_planes.clear();
        for (const auto& jp : cfg["planes"]) m_planes.push_back(jp.asInt());
    }
    m_tick0 = get(cfg, "tick0", m_tick0);
    m_nticks = get(cfg, "nticks", m_nticks);
    m_tick_span = get(cfg, "tick_span", m_tick_span);
    m_input_scale = get(cfg, "input_scale", m_input_scale);
    m_input_offset = get(cfg, "input_offset", m_input_offset);
    m_active_threshold = get(cfg, "active_threshold", m_active_threshold);
    if (cfg.isMember("view_norm")) {
        m_view_norm.clear();
        for (const auto& jv : cfg["view_norm"]) {
            m_view_norm.push_back({jv[0].asDouble(), jv[1].asDouble()});
        }
    }
    m_min_canvas = get(cfg, "min_canvas", m_min_canvas);
    m_bbox_pad = get(cfg, "bbox_pad", m_bbox_pad);
    m_feature_dim = get(cfg, "feature_dim", m_feature_dim);
    m_max_dense_pixels = get<Json::Int64>(cfg, "max_dense_pixels", m_max_dense_pixels);
    m_halo = get(cfg, "halo", m_halo);
    m_store_half = get(cfg, "store_half", m_store_half);
    if (cfg.isMember("provenance")) {
        m_provenance = cfg["provenance"];
    }

    if (m_input_tag.empty()) {
        log->critical("input_tag is empty");
        THROW(ValueError() << errmsg{"FMFeatureExtract: input_tag is empty"});
    }
    if (m_tick_span <= 0 or m_nticks <= 0 or m_nticks % m_tick_span != 0) {
        log->critical("nticks={} must be a positive multiple of tick_span={}", m_nticks, m_tick_span);
        THROW(ValueError() << errmsg{"FMFeatureExtract: nticks must be a positive multiple of tick_span"});
    }
    if (m_bbox_pad < 1) m_bbox_pad = 1;
    if (m_halo < 0) m_halo = 0;

    m_anode = Factory::find_tn<IAnodePlane>(m_anode_tn);
    m_forward = Factory::find_tn<ITensorForward>(m_forward_tn);

    m_rows.clear();
    for (int plane : m_planes) {
        if (plane < 0 or plane >= (int) m_view_norm.size()) {
            log->critical("plane {} has no view_norm entry ({} given)", plane, m_view_norm.size());
            THROW(ValueError() << errmsg{"FMFeatureExtract: plane without view_norm"});
        }
        std::set<int> idents;
        for (const auto& ich : Aux::plane_channels(m_anode, plane)) {
            idents.insert(ich->ident());
        }
        if (idents.empty()) {
            log->critical("anode={} plane={} has no channels", m_anode_tn, plane);
            THROW(ValueError() << errmsg{"FMFeatureExtract: plane without channels"});
        }
        PlaneRows pr;
        pr.plane = plane;
        pr.base = *idents.begin();
        pr.chlist.assign(idents.begin(), idents.end());
        pr.nrows = (int) pr.chlist.size();
        // The FM rebases each view by its first channel, which assumes contiguous idents.
        if (pr.chlist.back() - pr.base + 1 != pr.nrows) {
            log->critical("anode={} plane={} channel idents are not contiguous: [{},{}] for {} channels",
                          m_anode_tn, plane, pr.base, pr.chlist.back(), pr.nrows);
            THROW(ValueError() << errmsg{"FMFeatureExtract: non-contiguous channel idents"});
        }
        log->debug("anode={} plane={} rows={} channels=[{},{}] view_norm=({},{})",
                   m_anode_tn, plane, pr.nrows, pr.base, pr.chlist.back(),
                   m_view_norm[plane][0], m_view_norm[plane][1]);
        m_rows.push_back(pr);
    }
    log->debug("forward={} tag={} tick_span={} scale={} offset={} threshold={} min_canvas={} bbox_pad={} "
               "feature_dim={} max_dense_pixels={} halo={} store_half={}",
               m_forward_tn, m_input_tag, m_tick_span, m_input_scale, m_input_offset, m_active_threshold,
               m_min_canvas, m_bbox_pad, m_feature_dim, m_max_dense_pixels, m_halo, m_store_half);
}

namespace {
    int pad_to(int n, int mult)
    {
        if (mult <= 1) return n;
        return ((n + mult - 1) / mult) * mult;
    }
}

bool Pytorch::FMFeatureExtract::operator()(const input_pointer& in, output_pointer& out)
{
    out = nullptr;
    if (!in) {
        log->debug("EOS at call={}", m_count);
        return true;
    }
    const auto tstart = std::chrono::steady_clock::now();

    auto traces = Aux::tagged_traces(in, m_input_tag);
    if (traces.empty()) {
        log->warn("call={} frame={} has no traces tagged \"{}\"", m_count, in->ident(), m_input_tag);
    }

    const int nslices = m_nticks / m_tick_span;
    const int C = m_feature_dim;

    auto itv = std::make_shared<ITensor::vector>();
    Configuration md;
    md["frame_ident"] = in->ident();
    md["frame_time"] = in->time();
    md["frame_tick"] = in->tick();
    md["anode"] = m_anode->ident();
    md["tick0"] = m_tick0;
    md["nticks"] = m_nticks;
    md["tick_span"] = m_tick_span;
    md["input_tag"] = m_input_tag;
    md["input_scale"] = m_input_scale;
    md["input_offset"] = m_input_offset;
    md["active_threshold"] = m_active_threshold;
    md["face_layout"] = "channel";
    md["min_canvas"] = m_min_canvas;
    md["bbox_pad"] = m_bbox_pad;
    md["feature_dim"] = C;
    md["store_half"] = m_store_half;
    md["max_dense_pixels"] = (Json::Int64) m_max_dense_pixels;
    md["halo"] = m_halo;
    md["provenance"] = m_provenance;
    md["planes"] = Json::arrayValue;

    for (const auto& pr : m_rows) {
        const auto tplane = std::chrono::steady_clock::now();

        // 1. rows x ticks, charge summed per tick_span slice (double accumulation, then float).
        Array::array_xxf arr = Array::array_xxf::Zero(pr.nrows, m_nticks);
        {
            auto chlist = pr.chlist;
            Aux::fill(arr, traces, chlist.begin(), chlist.end(), m_tick0);
        }
        std::vector<float> q((size_t) pr.nrows * nslices, 0.0f);
        for (int r = 0; r < pr.nrows; ++r) {
            for (int s = 0; s < nslices; ++s) {
                double sum = 0.0;
                const int t0 = s * m_tick_span;
                for (int t = 0; t < m_tick_span; ++t) {
                    sum += arr(r, t0 + t);
                }
                q[(size_t) r * nslices + s] = (float) (sum * m_input_scale + m_input_offset);
            }
        }

        // 2. active pixels in (row, slice) order == (channel, slice) order, and their bounding box.
        std::vector<int> arow, acol;
        int r0 = std::numeric_limits<int>::max(), r1 = -1, k0 = std::numeric_limits<int>::max(), k1 = -1;
        for (int r = 0; r < pr.nrows; ++r) {
            for (int s = 0; s < nslices; ++s) {
                if (q[(size_t) r * nslices + s] > m_active_threshold) {
                    arow.push_back(r);
                    acol.push_back(s);
                    r0 = std::min(r0, r); r1 = std::max(r1, r);
                    k0 = std::min(k0, s); k1 = std::max(k1, s);
                }
            }
        }
        const size_t N = arow.size();

        Configuration pmd;
        pmd["plane"] = pr.plane;
        pmd["base_channel"] = pr.base;
        pmd["n_active"] = (Json::UInt64) N;
        pmd["view_norm"] = Json::arrayValue;
        pmd["view_norm"].append(m_view_norm[pr.plane][0]);
        pmd["view_norm"].append(m_view_norm[pr.plane][1]);

        std::vector<float> feat((size_t) N * C, 0.0f);
        int h = 0, w = 0, ntiles = 0;
        if (N) {
            // 3. the log map (float32 arithmetic, as torch applied it in training).
            const float m = (float) m_view_norm[pr.plane][0];
            const float y0 = (float) std::log10(m_view_norm[pr.plane][0]);
            const float sc = (float) (2.0 / (std::log10(m_view_norm[pr.plane][1] + m_view_norm[pr.plane][0]) - std::log10(m_view_norm[pr.plane][0])));
            std::vector<float> val(N);
            for (size_t i = 0; i < N; ++i) {
                const float x = q[(size_t) arow[i] * nslices + acol[i]];
                val[i] = (std::log10(x + m) - y0) * sc - 1.0f;
            }

            // 4. canvas geometry: tight bbox, floored, padded.
            h = pad_to(std::max(r1 - r0 + 1, m_min_canvas), m_bbox_pad);
            w = pad_to(std::max(k1 - k0 + 1, m_min_canvas), m_bbox_pad);

            // tiles along the slice axis: core [a,b) with halo on both sides.
            std::vector<std::pair<int, int>> cores;
            if (m_max_dense_pixels > 0 and (long) h * w > m_max_dense_pixels) {
                const int core = std::max(1, (int) (m_max_dense_pixels / h) - 2 * m_halo);
                for (int a = 0; a < w; a += core) {
                    cores.emplace_back(a, std::min(w, a + core));
                }
            }
            else {
                cores.emplace_back(0, w);
            }
            ntiles = (int) cores.size();

            // 5. forward per tile, gather the core columns.
            for (const auto& [a, b] : cores) {
                const int ea = std::max(0, a - m_halo);
                const int eb = std::min(w, b + m_halo);
                const int tw = pad_to(std::max(eb - ea, m_min_canvas), m_bbox_pad);
                std::vector<float> tile((size_t) 2 * h * tw, 0.0f);
                for (size_t i = 0; i < N; ++i) {
                    const int c = acol[i] - k0;
                    if (c < ea or c >= eb) continue;
                    const size_t idx = (size_t) (arow[i] - r0) * tw + (c - ea);
                    tile[idx] = val[i];
                    tile[(size_t) h * tw + idx] = 1.0f;
                }
                ITensor::shape_t ishape = {1, 2, (size_t) h, (size_t) tw};
                auto iten = std::make_shared<Aux::SimpleTensor>(ishape, tile.data());
                auto iset = std::make_shared<Aux::SimpleTensorSet>(
                    in->ident(), Configuration(), std::make_shared<ITensor::vector>(ITensor::vector{iten}));
                auto oset = m_forward->forward(iset);
                if (!oset or !oset->tensors() or oset->tensors()->empty()) {
                    log->critical("call={} plane={} forward \"{}\" returned nothing", m_count, pr.plane, m_forward_tn);
                    THROW(RuntimeError() << errmsg{"FMFeatureExtract: forward returned nothing"});
                }
                auto oten = oset->tensors()->front();
                const auto oshape = oten->shape();
                if (oshape.size() != 4 or oshape[0] != 1 or (int) oshape[1] != C or (int) oshape[2] != h or (int) oshape[3] != tw
                    or oten->dtype() != "f4") {
                    std::string got;
                    for (auto s : oshape) got += std::to_string(s) + ",";
                    log->critical("call={} plane={} forward reply shape/dtype mismatch: got [{}] {} expected [1,{},{},{}] f4",
                                  m_count, pr.plane, got, oten->dtype(), C, h, tw);
                    THROW(RuntimeError() << errmsg{"FMFeatureExtract: forward reply shape mismatch"});
                }
                const float* od = reinterpret_cast<const float*>(oten->data());
                for (size_t i = 0; i < N; ++i) {
                    const int c = acol[i] - k0;
                    if (c < a or c >= b) continue;
                    const size_t row = arow[i] - r0;
                    const size_t col = c - ea;
                    float* fi = feat.data() + i * C;
                    for (int ch = 0; ch < C; ++ch) {
                        fi[ch] = od[((size_t) ch * h + row) * tw + col];
                    }
                }
            }
            pmd["bbox"] = Json::arrayValue;
            pmd["bbox"].append(r0); pmd["bbox"].append(r1); pmd["bbox"].append(k0); pmd["bbox"].append(k1);
            pmd["canvas"] = Json::arrayValue;
            pmd["canvas"].append(h); pmd["canvas"].append(w);
            pmd["tiles"] = ntiles;
        }

        // 6. output tensors.
        std::vector<int32_t> coords(2 * N);
        for (size_t i = 0; i < N; ++i) {
            coords[2 * i] = pr.base + arow[i];
            coords[2 * i + 1] = acol[i];
        }
        Configuration cmd = pmd;
        cmd["name"] = "coords";
        ITensor::shape_t cshape = {N, 2};
        itv->push_back(std::make_shared<Aux::SimpleTensor>(cshape, coords.data(), cmd));

        Configuration fmd = pmd;
        ITensor::shape_t fshape = {N, (size_t) C};
        if (m_store_half) {
            std::vector<uint16_t> bits(feat.size());
            for (size_t i = 0; i < feat.size(); ++i) {
                bits[i] = c10::Half(feat[i]).x;
            }
            fmd["name"] = "feat_half";
            itv->push_back(std::make_shared<Aux::SimpleTensor>(fshape, bits.data(), fmd));
        }
        else {
            fmd["name"] = "feat";
            itv->push_back(std::make_shared<Aux::SimpleTensor>(fshape, feat.data(), fmd));
        }
        md["planes"].append(pmd);

        const double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - tplane).count();
        log->debug("call={} frame={} anode={} plane={} n_active={} bbox=[{},{},{},{}] canvas={}x{} tiles={} ms={:.0f}",
                   m_count, in->ident(), m_anode->ident(), pr.plane, N, r0, r1, k0, k1, h, w, ntiles, ms);
    }

    out = std::make_shared<Aux::SimpleTensorSet>(in->ident(), md, itv);
    const double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - tstart).count();
    log->debug("call={} frame={} anode={} planes={} tensors={} ms={:.0f}",
               m_count, in->ident(), m_anode->ident(), m_rows.size(), itv->size(), ms);
    ++m_count;
    return true;
}
