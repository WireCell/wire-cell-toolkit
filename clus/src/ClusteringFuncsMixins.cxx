
#include "WireCellClus/ClusteringFuncsMixins.h"
#include "WireCellUtil/NamedFactory.h"

using namespace WireCell;
using namespace WireCell::Clus;

void NeedDV:: configure(const WireCell::Configuration &cfg)
{
    auto tn = get<std::string>(cfg, "detector_volumes", "DetectorVolumes");
    m_dv = Factory::find_tn<IDetectorVolumes>(tn);
}
void NeedRecombModel:: configure(const WireCell::Configuration &cfg)
{
    auto tn = get<std::string>(cfg, "recombination_model", "BoxRecombination");
    m_recomb_model = Factory::find_tn<IRecombinationModel>(tn);
}

void NeedFiducial:: configure(const WireCell::Configuration &cfg)
{
    auto tn = get<std::string>(cfg, "fiducial", "DetectorVolumes");
    m_fiducial = Factory::find_tn<IFiducial>(tn);
}

void NeedPCTS::configure(const WireCell::Configuration &cfg)
{
    auto tn = get<std::string>(cfg, "pc_transforms", "PCTransformSet");
    m_pcts = Factory::find_tn<IPCTransformSet>(tn);
}

NeedScope::NeedScope(const std::string &pcname,
                     const std::vector<std::string> &coords,
                     size_t depth)
    : m_scope{pcname, coords, depth} 
{
}

void NeedScope::configure(const WireCell::Configuration &cfg) 
{
    m_scope.pcname = get(cfg, "pc_name", m_scope.pcname);
    m_scope.coords = get(cfg, "coords", m_scope.coords);
}
