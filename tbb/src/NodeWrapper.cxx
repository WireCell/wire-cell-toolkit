#include "WireCellTbb/NodeWrapper.h"
#include "WireCellUtil/Type.h"

// Off by default; DataFlowGraph turns it on when a timeline file is configured.
bool WireCellTbb::NodeInfo::s_collect_intervals = false;

std::ostream& WireCellTbb::operator<<(std::ostream& os, const WireCellTbb::NodeInfo& info)
{
    WireCell::INode::pointer n = info.inode();

    std::string tname = WireCell::type(*(n.get()));

    const double rtmax = std::chrono::duration_cast<std::chrono::milliseconds>(info.max_runtime()).count();
    const double rttot = std::chrono::duration_cast<std::chrono::milliseconds>(info.runtime()).count();
    const size_t num = info.calls();
    double rtmean = 0;

    if (num) rtmean = rttot/num;

    os 
       << "calls=" << num << " "
       << "time=" << rttot << " "
       << "mean=" << rtmean << " "
       << "max=" << rtmax << " [wall-ms] "
       << "core=" << info.coretime() << " [s] " 
       << "[" << tname << "] \"" << info.instance_name() << "\"";

    return os;
}
