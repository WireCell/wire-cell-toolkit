// Test the "any" interface of the four typed IHydraNode adapters.
//
// The IHydraNode contract says elements a body leaves on its input queues
// persist to the next call.  Callers of the "any" interface (eg the Pgraph
// Hydra wrapper) learn what the body consumed only from the any input queues,
// so each adapter must reflect the body's pops back into them.  Without that,
// Pgraph never pops a hydra's input and the node fires forever.

#include "WireCellIface/IHydraNode.h"
#include "WireCellIface/IDepo.h"
#include "WireCellUtil/doctest.h"

#include <boost/any.hpp>
#include <tuple>

using namespace WireCell;

namespace {

    struct TestDepo : public IDepo {
        int ident{0};
        Point point;
        TestDepo(int ident) : ident(ident) {}
        virtual ~TestDepo() {}
        virtual int id() const { return ident; }
        virtual int pdg() const { return 0; }
        virtual IDepo::pointer prior() const { return nullptr; }
        virtual const Point& pos() const { return point; }
        virtual double time() const { return 0.0; }
        virtual double charge() const { return 0.0; }
        virtual double energy() const { return 0.0; }
    };

    // Each body consumes only the front element of each input queue and
    // forwards it to output 0, leaving the rest for a later call.

    struct FrontTT : public IHydraNodeTT<std::tuple<IDepo, IDepo>, std::tuple<IDepo>> {
        virtual ~FrontTT() {}
        virtual std::string signature() { return typeid(*this).name(); }
        virtual bool operator()(input_queues& iqs, output_queues& oqs)
        {
            auto& q0 = std::get<0>(iqs);
            auto& q1 = std::get<1>(iqs);
            if (!q0.empty()) { std::get<0>(oqs).push_back(q0.front()); q0.pop_front(); }
            if (!q1.empty()) { std::get<0>(oqs).push_back(q1.front()); q1.pop_front(); }
            return true;
        }
    };

    struct FrontVV : public IHydraNodeVV<IDepo, IDepo, 2, 1> {
        virtual ~FrontVV() {}
        virtual std::string signature() { return typeid(*this).name(); }
        virtual bool operator()(input_queues& iqs, output_queues& oqs)
        {
            oqs.resize(1);
            for (auto& q : iqs) {
                if (!q.empty()) { oqs[0].push_back(q.front()); q.pop_front(); }
            }
            return true;
        }
    };

    struct FrontVT : public IHydraNodeVT<IDepo, std::tuple<IDepo>, 2> {
        virtual ~FrontVT() {}
        virtual std::string signature() { return typeid(*this).name(); }
        virtual bool operator()(input_queues& iqs, output_queues& oqs)
        {
            for (auto& q : iqs) {
                if (!q.empty()) { std::get<0>(oqs).push_back(q.front()); q.pop_front(); }
            }
            return true;
        }
    };

    struct FrontTV : public IHydraNodeTV<std::tuple<IDepo, IDepo>, IDepo, 1> {
        virtual ~FrontTV() {}
        virtual std::string signature() { return typeid(*this).name(); }
        virtual bool operator()(input_queues& iqs, output_queues& oqs)
        {
            oqs.resize(1);
            auto& q0 = std::get<0>(iqs);
            auto& q1 = std::get<1>(iqs);
            if (!q0.empty()) { oqs[0].push_back(q0.front()); q0.pop_front(); }
            if (!q1.empty()) { oqs[0].push_back(q1.front()); q1.pop_front(); }
            return true;
        }
    };

    int any_id(const boost::any& a)
    {
        return boost::any_cast<IDepo::pointer>(a)->id();
    }

    // Two inputs with two depos each: {0,1} and {10,11}.
    IHydraNodeBase::any_queue_vector make_inputs()
    {
        IHydraNodeBase::any_queue_vector inqs(2);
        for (int port = 0; port < 2; ++port) {
            for (int ind = 0; ind < 2; ++ind) {
                IDepo::pointer depo = std::make_shared<TestDepo>(10 * port + ind);
                inqs[port].push_back(depo);
            }
        }
        return inqs;
    }

    void check_front_consumed(IHydraNodeBase& node)
    {
        auto inqs = make_inputs();
        IHydraNodeBase::any_queue_vector outqs;
        REQUIRE(node(inqs, outqs));

        // The body's pops are reflected in the any input queues.
        REQUIRE(inqs.size() == 2);
        REQUIRE(inqs[0].size() == 1);
        REQUIRE(inqs[1].size() == 1);
        CHECK(any_id(inqs[0].front()) == 1);
        CHECK(any_id(inqs[1].front()) == 11);

        REQUIRE(outqs.size() == 1);
        REQUIRE(outqs[0].size() == 2);
        CHECK(any_id(outqs[0][0]) == 0);
        CHECK(any_id(outqs[0][1]) == 10);
    }
}

TEST_CASE("iface hydra any interface reflects consumed input")
{
    SUBCASE("TT") { FrontTT node; check_front_consumed(node); }
    SUBCASE("VV") { FrontVV node; check_front_consumed(node); }
    SUBCASE("VT") { FrontVT node; check_front_consumed(node); }
    SUBCASE("TV") { FrontTV node; check_front_consumed(node); }
}
