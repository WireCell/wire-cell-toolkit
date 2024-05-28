#include "WireCellGen/DeposOrBust.h"
#include "WireCellIface/IFrame.h"

#include "WireCellUtil/NamedFactory.h"

WIRECELL_FACTORY(DeposOrBust,
                 WireCell::Gen::DeposOrBust,
                 WireCell::IDepos2DeposOrFrame)


using namespace WireCell;

Gen::DeposOrBust::DeposOrBust()
    : Aux::Logger("DeposOrBust", "glue")
    { }
Gen::DeposOrBust::~DeposOrBust() { }

// When I get an empty deposet, this is what comes out of port 1.
class EmptyFrame : public IFrame {
    int _ident{0};
  public:
    EmptyFrame(int ident) : _ident(ident) {}
    ~EmptyFrame() {}

    virtual const tag_list_t& frame_tags() const {
        static tag_list_t dummy;
        return dummy;
    }

    virtual const tag_list_t& trace_tags() const {
        static tag_list_t dummy;
        return dummy;
    }
    virtual const trace_list_t& tagged_traces(const tag_t&) const {
        static trace_list_t dummy;
        return dummy;
    }
    virtual const trace_summary_t& trace_summary(const tag_t&) const {
        static trace_summary_t dummy;
        return dummy;
    }
    virtual ITrace::shared_vector traces() const {
        return nullptr;
    }

    virtual int ident() const { return _ident; }
    virtual double time() const { return 0.0; }
    virtual double tick() const { return 0.0; }

};


bool Gen::DeposOrBust::operator()(input_queues& iqs, output_queues& oqs)
{
    auto& iq = std::get<0>(iqs);

    while (! iq.empty()) {
        auto ds = iq.front();
        iq.pop_front();
        
        if (!ds) {              // EOS on every ouput
            log->debug("EOS sync");
            get<0>(oqs).push_back(nullptr);
            get<1>(oqs).push_back(nullptr);
            continue;
        }

        if (ds->depos()->empty()) {
            log->debug("DeposOrBust: got empty deposet");
            get<1>(oqs).push_back(std::make_shared<EmptyFrame>(ds->ident()));
        }
        else {
            log->debug("DeposOrBust: forward deposet");
            get<0>(oqs).push_back(ds);
        }
    }
    return true;
}


