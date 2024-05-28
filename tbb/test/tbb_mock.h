#ifndef WIRECELLTBB_TBB_MOCK
#define WIRECELLTBB_TBB_MOCK

#include "WireCellIface/IDepoSource.h"
#include "WireCellIface/IDrifter.h"
#include "WireCellIface/IDepoSink.h"
#include "WireCellAux/Logger.h"
#include "WireCellAux/SimpleDepo.h"
#include "WireCellUtil/Units.h"

#include <iostream>
#include <chrono>
#include <thread>

namespace WireCellTbb {

    using namespace std::chrono_literals;

    class MockDepoSource : public WireCell::IDepoSource, public WireCell::Aux::Logger {
        int m_count;
        const int m_maxdepos;

      public:
        MockDepoSource(int maxdepos = 1000)
            : WireCell::Aux::Logger("MockDepoSource","test")
            , m_count(0)
            , m_maxdepos(maxdepos)
        {
        }
        virtual ~MockDepoSource() {}

        virtual bool operator()(output_pointer& out)
        {
            out = nullptr;

            if (m_count == m_maxdepos) {
                std::cerr << "Source: EOS\n";
                ++m_count;
                return true;
            }
            if (m_count > m_maxdepos) {
                return false;
            }
            std::this_thread::sleep_for(1ms);
            double dist = m_count * WireCell::units::millimeter;
            double time = m_count * WireCell::units::microsecond;
            WireCell::Point pos(dist, dist, dist);
            out = WireCell::IDepo::pointer(new WireCell::Aux::SimpleDepo(time, pos));
            std::cerr << "Source: " << out->time() / WireCell::units::millimeter << std::endl;
            return true;
        }
    };

    class MockDrifter : public WireCell::IDrifter, public WireCell::Aux::Logger {
        std::deque<input_pointer> m_depos;

      public:
        MockDrifter ()
            : WireCell::Aux::Logger("MockDrifter","test")
        {
        }
        virtual ~MockDrifter() {}

        virtual bool operator()(const input_pointer& in, output_queue& outq)
        {
            std::this_thread::sleep_for(2ms);

            m_depos.push_back(in);

            // simulate some buffering condition
            size_t n_to_keep = 2;
            if (!in) {
                n_to_keep = 0;
            }

            while (m_depos.size() > n_to_keep) {
                auto depo = m_depos.front();
                m_depos.pop_front();
                outq.push_back(depo);
                if (!depo) {
                    std::cerr << "Drift: EOS\n";
                }
                else {
                    std::cerr << "Drift: " << depo->time() / WireCell::units::microsecond << std::endl;
                }
            }

            return true;
        }
    };

    class MockDepoSink : public WireCell::IDepoSink, public WireCell::Aux::Logger {
      public:
        MockDepoSink()
            : WireCell::Aux::Logger("MockDepoSink","test")
        {
        }
        virtual ~MockDepoSink() {}
        virtual bool operator()(const input_pointer& depo)
        {
            //std::this_thread::sleep_for(1ms);

            if (!depo) {
                std::cerr << "Sink: EOS\n";
            }
            else {
                std::cerr << "Sink: " << depo->time() / WireCell::units::microsecond << std::endl;
            }
            return true;
        }
    };

}  // namespace WireCellTbb
#endif
