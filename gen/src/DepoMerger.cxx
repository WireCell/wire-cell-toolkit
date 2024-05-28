#include "WireCellGen/DepoMerger.h"

#include "WireCellUtil/NamedFactory.h"

WIRECELL_FACTORY(DepoMerger, WireCell::Gen::DepoMerger, WireCell::IDepoMerger, WireCell::IConfigurable)

using namespace WireCell;

Gen::DepoMerger::DepoMerger()
  : m_nin0(0)
  , m_nin1(0)
  , m_nout(0)
  , m_eos(false)
{
}
Gen::DepoMerger::~DepoMerger() {}

bool Gen::DepoMerger::operator()(input_queues& inqs, output_queues& outqs)
{
    if (m_eos) {  // already closed off all our outputs
        return false;
    }

    auto& inq0 = get<0>(inqs);
    auto& inq1 = get<1>(inqs);

    if (inq0.empty() or inq1.empty()) {
        std::cerr << "DepoMerger: called empty input\n";
        return false;
    }

    auto& outq = get<0>(outqs);

    // std::cerr << "DepoMerger queue sizes on input: in:["<< inq0.size()<<","<<inq1.size()<<"], "
    //          << "out:["<<outq.size()<<"]\n";

    IDepo::pointer d0 = inq0.front();
    IDepo::pointer d1 = inq1.front();

    if (d0 and d1) {
        double t0 = d0->time();  // keep the newest one
        double t1 = d1->time();  // which may be both if they coincide
        if (t0 <= t1) {
            ++m_nout;
            ++m_nin0;
            outq.push_back(d0);
            inq0.pop_front();
            // std::cerr << "DepoMerger: stream 0 output: t0="<<t0<<", t1="<<t1<<"\n";
        }
        if (t1 <= t0) {
            ++m_nout;
            ++m_nin1;
            outq.push_back(d1);
            inq1.pop_front();
            // std::cerr << "DepoMerger: stream 1 output: t0="<<t0<<", t1="<<t1<<"\n";
        }
        return true;
    }

    if (d0) {  // d1 is eos
        ++m_nout;
        ++m_nin0;
        outq.push_back(d0);
        inq0.pop_front();
        // std::cerr << "DepoMerger: stream 0 only output: nout="<<m_nout<<"\n";
        return true;
    }

    if (d1) {  // d0 is eos
        ++m_nout;
        ++m_nin1;
        outq.push_back(d1);
        inq1.pop_front();
        // std::cerr << "DepoMerger: stream 1 only output: nout="<<m_nout<<"\n";
        return true;
    }

    // both are eos for the first time.
    m_eos = true;
    outq.push_back(nullptr);
    std::cerr << "DepoMerger: global EOS: in: " << m_nin0 << " + " << m_nin1 << ", out: " << m_nout << " depos\n";

    return true;
}

void Gen::DepoMerger::configure(const WireCell::Configuration& cfg) {}

WireCell::Configuration Gen::DepoMerger::default_configuration() const
{
    Configuration cfg;
    return cfg;
}
