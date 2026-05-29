#include "WireCellMatch/OpflashToFlashPCs.h"

#include "WireCellAux/TensorDMcommon.h"    // as_tensorset
#include "WireCellAux/TensorDMpointtree.h" // as_pctree, as_tensors

#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/PointCloudArray.h"
#include "WireCellUtil/PointTree.h"
#include "WireCellUtil/String.h"

WIRECELL_FACTORY(OpflashToFlashPCs,
                 WireCell::Match::OpflashToFlashPCs,
                 WireCell::INamed,
                 WireCell::ITensorSetFanin,
                 WireCell::IConfigurable)

using namespace WireCell;

Match::OpflashToFlashPCs::OpflashToFlashPCs()
    : Aux::Logger("OpflashToFlashPCs", "match")
{
}

Match::OpflashToFlashPCs::~OpflashToFlashPCs() = default;

void Match::OpflashToFlashPCs::configure(const WireCell::Configuration& cfg)
{
    m_inpath = get(cfg, "inpath", m_inpath);
    m_nchan  = get(cfg, "nchan", m_nchan);
    log->debug("OpflashToFlashPCs: inpath={} nchan={}", m_inpath, m_nchan);
}

WireCell::Configuration Match::OpflashToFlashPCs::default_configuration() const
{
    Configuration cfg;
    cfg["inpath"] = m_inpath;
    cfg["nchan"]  = m_nchan;
    return cfg;
}

std::vector<std::string> Match::OpflashToFlashPCs::input_types()
{
    const std::string tname = std::string(typeid(input_type).name());
    return std::vector<std::string>(m_multiplicity, tname);
}

bool Match::OpflashToFlashPCs::operator()(const input_vector& invec, output_pointer& out)
{
    out = nullptr;

    // Inputs are framework-synced: any EOS is EOS.
    for (const auto& in : invec) {
        if (!in) {
            log->debug("EOS at call {}", m_count++);
            return true;
        }
    }

    const auto& tree_ts = invec[0]; // pctree (live + dead)
    const auto& data_ts = invec[1]; // opflash matrix tensor set

    const int ident = tree_ts->ident();
    std::string inpath = m_inpath;
    if (inpath.find("%") != std::string::npos) inpath = String::format(inpath, ident);

    const auto& tree_tens = *tree_ts->tensors();
    auto root_live = Aux::TensorDM::as_pctree(tree_tens, inpath + "/live");
    if (!root_live) {
        raise<ValueError>("OpflashToFlashPCs: no live pctree at '%s'", inpath);
    }

    // Canonical optical PCs (same schema as root/UbooneClusterSource).
    std::vector<double> ftime, ftmin, ftmax, fval;
    std::vector<int>    fident, ftype;
    std::vector<int>    lid;
    std::vector<double> lt, lq, lerr;
    std::vector<int>    fl_flash, fl_light;

    const auto& data_tens = data_ts->tensors();
    if (data_tens && !data_tens->empty()) {
        const auto& ten = data_tens->at(0);
        if (ten->shape().size() != 2) {
            raise<ValueError>("OpflashToFlashPCs: tensor dim %d != 2", ten->shape().size());
        }
        const size_t nflash = ten->shape()[0];
        const size_t ncol   = ten->shape()[1];
        if ((int)ncol != m_nchan + 1) {
            raise<ValueError>("OpflashToFlashPCs: ncol %d != nchan+1 (%d)", (int)ncol, m_nchan + 1);
        }
        const double* M = (const double*) ten->data(); // row-major [nflash][ncol]

        for (size_t r = 0; r < nflash; ++r) {
            const double time = M[r * ncol + 0];
            double sum = 0;
            for (int c = 0; c < m_nchan; ++c) {
                const double pe = M[r * ncol + 1 + c];
                sum += pe;
                if (pe != 0.0) {
                    // flashlight join (append before light so light index = lid.size()).
                    fl_flash.push_back((int)r);
                    fl_light.push_back((int)lid.size());
                    lid.push_back(c);
                    lt.push_back(time);
                    lq.push_back(pe);
                    lerr.push_back(0.0);
                }
            }
            ftime.push_back(time);
            ftmin.push_back(time); // SBND matrix has no per-flash window
            ftmax.push_back(time);
            fval.push_back(sum);
            fident.push_back((int)r);
            ftype.push_back(0);
        }
        log->debug("OpflashToFlashPCs: ident {} {} flashes, {} light entries",
                   ident, nflash, lid.size());
    }
    else {
        log->warn("OpflashToFlashPCs: ident {} no tensor on port 1; attaching empty flash PCs", ident);
    }

    PointCloud::Dataset flash_ds;
    flash_ds.add("time",  PointCloud::Array(ftime));
    flash_ds.add("tmin",  PointCloud::Array(ftmin));
    flash_ds.add("tmax",  PointCloud::Array(ftmax));
    flash_ds.add("value", PointCloud::Array(fval));
    flash_ds.add("ident", PointCloud::Array(fident));
    flash_ds.add("type",  PointCloud::Array(ftype));

    PointCloud::Dataset light_ds;
    light_ds.add("ident", PointCloud::Array(lid));
    light_ds.add("time",  PointCloud::Array(lt));
    light_ds.add("value", PointCloud::Array(lq));
    light_ds.add("error", PointCloud::Array(lerr));

    PointCloud::Dataset flashlight_ds;
    flashlight_ds.add("flash", PointCloud::Array(fl_flash));
    flashlight_ds.add("light", PointCloud::Array(fl_light));

    auto& lpcs = root_live->value.local_pcs();
    lpcs["flash"]      = std::move(flash_ds);
    lpcs["light"]      = std::move(light_ds);
    lpcs["flashlight"] = std::move(flashlight_ds);

    // Re-serialize live (now carrying the optical PCs) and pass /dead through.
    ITensor::vector outtens;
    auto tens_live = Aux::TensorDM::as_tensors(*root_live, inpath + "/live");
    outtens.insert(outtens.end(), tens_live.begin(), tens_live.end());

    auto root_dead = Aux::TensorDM::as_pctree(tree_tens, inpath + "/dead");
    auto tens_dead = Aux::TensorDM::as_tensors(*root_dead, inpath + "/dead");
    outtens.insert(outtens.end(), tens_dead.begin(), tens_dead.end());

    out = Aux::TensorDM::as_tensorset(outtens, ident);
    ++m_count;
    return true;
}
