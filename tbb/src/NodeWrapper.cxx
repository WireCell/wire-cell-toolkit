#include "WireCellTbb/NodeWrapper.h"
#include "WireCellUtil/Type.h"

std::ostream& WireCellTbb::operator<<(std::ostream& os, const WireCellTbb::NodeInfo& info)
{
    WireCell::INode::pointer n = info.inode();
    auto tname = WireCell::demangle(n->signature());

    const double rtmax = std::chrono::duration_cast<std::chrono::milliseconds>(info.max_runtime()).count();
    const double rttot = std::chrono::duration_cast<std::chrono::milliseconds>(info.runtime()).count();
    const size_t num = info.calls();
    double rtmean = 0;
    if (num) rtmean = rttot/num;

    os <<  tname << " \"" << info.instance_name() << "\" "
       << "calls=" << num << " "
       << "time=" << rttot << " "
       << "mean=" << rtmean << " "
       << "max=" << rtmax << " [ms]";

    return os;
}
