#include "WireCellTruth/BlobDepoLabelTensor.h"
#include "WireCellTruth/BlobDepoLabel.h"
#include "WireCellTruth/Schema.h"

#include "WireCellAux/SimpleTensor.h"
#include "WireCellAux/SimpleTensorSet.h"
#include "WireCellAux/TensorDMcluster.h"

#include "WireCellUtil/GraphTools.h"
#include "WireCellUtil/String.h"
#include "WireCellUtil/NamedFactory.h"

WIRECELL_FACTORY(BlobDepoLabelTensor, WireCell::Truth::BlobDepoLabelTensor,
                 WireCell::INamed,
                 WireCell::IClusterDepoTensor, WireCell::IConfigurable)

using namespace WireCell;
using WireCell::GraphTools::mir;

Truth::BlobDepoLabelTensor::BlobDepoLabelTensor()
    : Aux::Logger("BlobDepoLabelTensor", "truth")
{
}
Truth::BlobDepoLabelTensor::~BlobDepoLabelTensor()
{
}

Configuration Truth::BlobDepoLabelTensor::default_configuration() const
{
    auto cfg = assoc_default_configuration(m_params);
    cfg["datapath"] = m_datapath;
    cfg["with_cluster"] = m_with_cluster;
    cfg["min_frac"] = m_min_frac;
    return cfg;
}

void Truth::BlobDepoLabelTensor::configure(const Configuration& cfg)
{
    assoc_configure(m_params, cfg);
    m_datapath = get(cfg, "datapath", m_datapath);
    m_with_cluster = get(cfg, "with_cluster", m_with_cluster);
    m_min_frac = get(cfg, "min_frac", m_min_frac);
}

static Configuration columns_md(const std::vector<std::string>& cols)
{
    Configuration jcols = Json::arrayValue;
    for (const auto& c : cols) {
        jcols.append(c);
    }
    return jcols;
}


ITensor::vector Truth::blobtruth_tensors(const ICluster::pointer& icluster,
                                         const DepoBlobAssocResult& res,
                                         const std::string& datapath,
                                         double min_frac, size_t* nrows_out)
{
    // Rows in graph vertex order; tracks in key order.
    std::vector<double> table;   // M x 6
    std::vector<double> summary; // B x 7
    size_t nrows = 0, nblobs = 0;
    const auto& gr = icluster->graph();
    for (const auto& vdesc : mir(boost::vertices(gr))) {
        if (gr[vdesc].code() != 'b') continue;
        const auto iblob = std::get<IBlob::pointer>(gr[vdesc].ptr);
        const auto& e = res.blobs.at(vdesc);
        const double bident = iblob->ident();
        const double face = iblob->face() ? iblob->face()->ident() : -1;
        const double sident = iblob->slice() ? iblob->slice()->ident() : -1;

        double q_top = 0.0;
        int tid_top = Truth::tid_ghost;
        size_t ntracks = 0;
        for (const auto& [tid, qe] : e.by_track) {
            const double q = qe.first;
            if (q <= 0.0) continue;
            if (min_frac > 0.0 && e.q_total > 0.0 && q/e.q_total < min_frac) continue;
            table.insert(table.end(), {bident, face, sident, (double)tid, q, qe.second});
            ++nrows;
            ++ntracks;
            if (q > q_top) { q_top = q; tid_top = tid; }
        }
        summary.insert(summary.end(), {bident, face, sident, e.q_total, q_top, (double)tid_top, (double)ntracks});
        ++nblobs;
    }

    ITensor::vector out;
    {
        Configuration md;
        md["datatype"] = Truth::blobtruth_datatype;
        md["datapath"] = datapath + "/blobtruth";
        md["ident"] = icluster->ident();
        md["columns"] = columns_md(Truth::blobtruth_columns());
        ITensor::shape_t shape{nrows, Truth::blobtruth_columns().size()};
        out.push_back(std::make_shared<Aux::SimpleTensor>(shape, table.data(), md));
    }
    {
        Configuration md;
        md["datatype"] = "blobtruth_sum";
        md["datapath"] = datapath + "/blobtruth_sum";
        md["ident"] = icluster->ident();
        md["columns"] = columns_md(Truth::blobtruth_sum_columns());
        md["ndepos"] = (Json::UInt64)res.ndepos;      // depos in the faces' sensitive volumes
        md["q_depos"] = res.q_depos;                   // their |charge| sum
        md["ndepos_window"] = (Json::UInt64)res.ndepos_window;   // of those, inside the slices' time span
        md["q_depos_window"] = res.q_depos_window;     // the completeness denominator
        md["nblobs"] = (Json::UInt64)nblobs;
        ITensor::shape_t shape{nblobs, Truth::blobtruth_sum_columns().size()};
        out.push_back(std::make_shared<Aux::SimpleTensor>(shape, summary.data(), md));
    }
    if (nrows_out) *nrows_out = nrows;
    return out;
}

bool Truth::BlobDepoLabelTensor::operator()(const input_tuple_type& intup,
                                            output_pointer& out)
{
    out = nullptr;
    auto icluster = std::get<0>(intup);
    auto ideposet = std::get<1>(intup);
    if (!icluster or !ideposet) {
        if (icluster) {
            log->warn("join node got non-EOS ICluster {} on EOS", icluster->ident());
        }
        if (ideposet) {
            log->warn("join node got non-EOS IDepoSet {} on EOS", ideposet->ident());
        }
        log->debug("EOS at call={}", m_count);
        ++m_count;
        return true;
    }

    auto res = depo_blob_assoc(icluster, *ideposet->depos(), m_params, log);
    const int ident = icluster->ident();
    const std::string datapath = String::format(m_datapath, ident);

    auto tv = std::make_shared<ITensor::vector>();
    if (m_with_cluster) {
        auto filled = filled_cluster(icluster, res);
        auto ctens = Aux::TensorDM::as_tensors(filled, datapath);
        tv->insert(tv->end(), ctens.begin(), ctens.end());
    }
    size_t nrows = 0;
    auto btens = blobtruth_tensors(icluster, res, datapath, m_min_frac, &nrows);
    tv->insert(tv->end(), btens.begin(), btens.end());
    const size_t nblobs = res.nblobs;

    log->debug("call={} cluster {} deposet {}: nblobs={} ndepos={} rows={}",
               m_count, ident, ideposet->ident(), nblobs, res.ndepos, nrows);

    Configuration setmd;
    setmd["ident"] = ident;
    out = std::make_shared<Aux::SimpleTensorSet>(ident, setmd, tv);
    ++m_count;
    return true;
}
