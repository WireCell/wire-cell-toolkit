#include "WireCellSpng/MaskRebaseliner.h"
#include "WireCellSpng/SimpleTorchTensor.h"
#include "WireCellSpng/TdmTools.h"
#include "WireCellSpng/Util.h"

#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/NamedFactory.h"

WIRECELL_FACTORY(SPNGMaskRebaseliner,
                 WireCell::SPNG::MaskRebaseliner,
                 WireCell::SPNG::MaskRebaseliner::fan_type,
                 WireCell::IConfigurable,
                 WireCell::INamed);

namespace WireCell::SPNG {

    MaskRebaseliner::MaskRebaseliner()
        : FaninBase<ITorchTensor>("MaskRebaseliner")
    {
    }

    void MaskRebaseliner::fanin_combine(const input_vector& inv, output_pointer& out)
    {
        if (inv.size() != 2) {
            raise<ValueError>("MaskRebaseliner: want 2 inputs (waveform, mask), got %d", inv.size());
        }
        auto wave = inv[0]->tensor();
        auto mask = inv[1]->tensor();
        if (wave.sizes() != mask.sizes()) {
            raise<ValueError>("MaskRebaseliner: waveform %s and mask %s differ in shape",
                              to_string(wave), to_string(mask));
        }
        log->debug("waveform {} mask {} enabled={} pad={} width={} keep_pad={} rows={}",
                   to_string(wave), to_string(mask), m_config.enabled,
                   m_config.anchor_pad, m_config.anchor_width, m_config.keep_pad, m_rows.size());

        torch::Tensor res;
        if (m_config.enabled) {
            res = rebaseline_masked_rows(wave, mask, m_rows, m_opts);
        }
        else {
            auto m = mask.scalar_type() == torch::kBool ? mask : mask.ne(0);
            res = torch::where(m, wave, torch::zeros_like(wave));
        }

        // Derive from the waveform only: the mask carries its own tag/datapath
        // which must not leak into the signal tensor.
        auto md = TDM::derive_metadata(Configuration(), inv[0]->metadata(),
                                       m_config.datapath_format, m_config.tag);
        out = std::make_shared<SimpleTorchTensor>(res, md);
        logit(out, "maskrebaseline");
    }

    void MaskRebaseliner::configure(const WireCell::Configuration& config)
    {
        FaninBase<ITorchTensor>::configure(config);
        HanaJsonCPP::from_json(m_config, config);

        m_opts.dim = m_config.dim;
        m_opts.anchor_pad = m_config.anchor_pad;
        m_opts.anchor_width = m_config.anchor_width;
        m_opts.min_roi_size = m_config.min_roi_size;
        m_opts.keep_pad = m_config.keep_pad;
        m_opts.remove_negative = m_config.remove_negative;

        m_rows.clear();
        for (const auto& r : m_config.rows) {
            if (r.size() != 2 || r[0] < 0 || r[1] < r[0]) {
                raise<ValueError>("MaskRebaseliner: bad row range, want [begin,end) with 0<=begin<=end");
            }
            m_rows.emplace_back(r[0], r[1]);
            log->debug("rebaseline rows [{},{})", r[0], r[1]);
        }
    }

    WireCell::Configuration MaskRebaseliner::default_configuration() const
    {
        auto cfg = FaninBase<ITorchTensor>::default_configuration();
        update(cfg, HanaJsonCPP::to_json(m_config));
        return cfg;
    }

}
