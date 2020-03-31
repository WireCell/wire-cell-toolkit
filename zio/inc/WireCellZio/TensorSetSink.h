/**
 * WireCell::ITensorSet -> zio::Message
 */

#ifndef WIRECELLZIO_ZIOTENSORSETSINK
#define WIRECELLZIO_ZIOTENSORSETSINK

#include "WireCellZio/FlowConfigurable.h"
#include "WireCellIface/ITensorSetSink.h"
#include "WireCellUtil/Logging.h"

namespace WireCell
{
namespace Zio
{

class TensorSetSink : public ITensorSetSink, public FlowConfigurable
{
public:
    TensorSetSink();
    virtual ~TensorSetSink();

    virtual bool operator()(const ITensorSet::pointer &in);

private:
    Log::logptr_t l;
    bool m_had_eos;

    virtual void post_configure();

};
} // namespace Zio
} // namespace WireCell
#endif // WIRECELLZIO_ZIOTENSORSETSINK