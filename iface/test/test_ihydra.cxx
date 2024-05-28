// this only tests the typed interface.

#include "WireCellIface/IHydraNode.h"
#include "WireCellIface/IDepoSource.h"
#include "WireCellUtil/Type.h"
#include <iostream>
#include <tuple>

using namespace WireCell;

struct DummyDepo : public IDepo {
    int ident{0};
    Point point;

    DummyDepo(int ident) : ident(ident) {}
    virtual ~DummyDepo() {}
    virtual int id() const { return ident; }
    virtual int pdg() const { return 0; }
    virtual IDepo::pointer prior() const { return nullptr; }
    virtual const Point& pos() const { return point; }
    virtual double time() const { return 0.0; }
    virtual double charge() const { return 0.0; }
    virtual double energy() const { return 0.0; }

};

class MyDepoSrc : public IDepoSource
{
    size_t count{0};

  public:
    virtual ~MyDepoSrc() {}

    virtual bool operator()(WireCell::IDepo::pointer& depo)
    {
        if (!count) {
            std::cerr << "Running instance of " << demangle(signature()) << "\n";

            std::cerr << "input types:\n";
            for (auto tn : input_types()) {
                std::cerr << "\t" << demangle(tn) << "\n";
            }
            
            std::cerr << "output types:\n";
            for (auto tn : output_types()) {
                std::cerr << "\t" << demangle(tn) << "\n";
            }
        }
        
        depo = std::make_shared<DummyDepo>(count++);

        return true;
    }
};

class MyHydraTT : public IHydraNodeTT<std::tuple<IDepo, IDepo>, std::tuple<IDepo>>
{
    bool first{true};

  public:

    using MyHydra = IHydraNodeTT<std::tuple<IDepo, IDepo>, std::tuple<IDepo>>;
    using MyHydra::input_queues;
    using MyHydra::output_queues;


    virtual ~MyHydraTT() {}

    virtual bool operator()(input_queues& iqs, output_queues& oqs)
    {
        if (first) {
            first = false;

            std::cerr << "Running instance of " << demangle(signature()) << "\n";
            std::cerr << "input_queues:\n\t" << demangle(typeid(input_queues).name()) << "\n";
            std::cerr << "output_queues:\n\t" << demangle(typeid(output_queues).name()) << "\n";

            std::cerr << "input types:\n";
            for (auto tn : input_types()) {
                std::cerr << "\t" << demangle(tn) << "\n";
            }

            std::cerr << "output types:\n";
            for (auto tn : output_types()) {
                std::cerr << "\t" << demangle(tn) << "\n";
            }
        }
        
        for (auto depo : std::get<0>(iqs)) {
            std::get<0>(oqs).push_back(depo);
        }
        for (auto depo : std::get<1>(iqs)) {
            std::get<0>(oqs).push_back(depo);
        }

        return true;
    }

    virtual std::string signature() { return typeid(MyHydraTT).name(); }
};

class MyHydraVV : public IHydraNodeVV<IDepo, IDepo, 2, 1>
{
    bool first{true};

  public:

    using MyHydra = IHydraNodeVV<IDepo, IDepo, 2, 1>;
    using MyHydra::input_queues;
    using MyHydra::output_queues;

    virtual ~MyHydraVV() {}

    virtual bool operator()(input_queues& iqs, output_queues& oqs)
    {
        if (first) {
            first = false;

            std::cerr << "Running instance of " << demangle(signature()) << "\n";
            std::cerr << "input_queues:\n\t" << demangle(typeid(input_queues).name()) << "\n";
            std::cerr << "output_queues:\n\t" << demangle(typeid(output_queues).name()) << "\n";

            std::cerr << "input types:\n";
            for (auto tn : input_types()) {
                std::cerr << "\t" << demangle(tn) << "\n";
            }

            std::cerr << "output types:\n";
            for (auto tn : output_types()) {
                std::cerr << "\t" << demangle(tn) << "\n";
            }
            assert(input_types().size() == 2);
            assert(output_types().size() == 1);
        }
        
        for (auto& depos : iqs) {
            for (auto& depo : depos) {
                oqs[0].push_back(depo);
            }
        }

        return true;
    }

    virtual std::string signature() { return typeid(MyHydraVV).name(); }
};


class MyHydraTV : public IHydraNodeTV<std::tuple<IDepo,IDepo>, IDepo, 1>
{
    bool first{true};

  public:

    using MyHydra = IHydraNodeTV<std::tuple<IDepo,IDepo>, IDepo, 1>;
    using MyHydra::input_queues;
    using MyHydra::output_queues;

    virtual ~MyHydraTV() {}

    virtual bool operator()(input_queues& iqs, output_queues& oqs)
    {
        if (first) {
            first = false;

            std::cerr << "Running instance of " << demangle(signature()) << "\n";
            std::cerr << "input_queues:\n\t" << demangle(typeid(input_queues).name()) << "\n";
            std::cerr << "output_queues:\n\t" << demangle(typeid(output_queues).name()) << "\n";

            std::cerr << "input types:\n";
            for (auto tn : input_types()) {
                std::cerr << "\t" << demangle(tn) << "\n";
            }

            std::cerr << "output types:\n";
            for (auto tn : output_types()) {
                std::cerr << "\t" << demangle(tn) << "\n";
            }
            assert(output_types().size() == 1);
        }
        
        for (auto depo : std::get<0>(iqs)) {
            oqs[0].push_back(depo);
        }
        for (auto depo : std::get<1>(iqs)) {
            oqs[0].push_back(depo);
        }

        return true;
    }

    virtual std::string signature() { return typeid(MyHydraTV).name(); }
};

class MyHydraVT : public IHydraNodeVT<IDepo, std::tuple<IDepo>, 2>
{
    bool first{true};

  public:

    using MyHydra = IHydraNodeVT<IDepo, std::tuple<IDepo>, 2>;
    using MyHydra::input_queues;
    using MyHydra::output_queues;

    virtual ~MyHydraVT() {}

    virtual bool operator()(input_queues& iqs, output_queues& oqs)
    {
        if (first) {
            first = false;

            std::cerr << "Running instance of " << demangle(signature()) << "\n";
            std::cerr << "input_queues:\n\t" << demangle(typeid(input_queues).name()) << "\n";
            std::cerr << "output_queues:\n\t" << demangle(typeid(output_queues).name()) << "\n";

            std::cerr << "input types:\n";
            for (auto tn : input_types()) {
                std::cerr << "\t" << demangle(tn) << "\n";
            }

            std::cerr << "output types:\n";
            for (auto tn : output_types()) {
                std::cerr << "\t" << demangle(tn) << "\n";
            }
            assert(input_types().size() == 2);
        }
        
        for (auto& depos : iqs) {
            for (auto& depo : depos) {
                get<0>(oqs).push_back(depo);
            }
        }

        return true;
    }

    virtual std::string signature() { return typeid(MyHydraVV).name(); }
};

void dump_queue(const std::deque<IDepo::pointer>& depos)
{
    std::cerr << "[" << depos.size() << "]";
    for (const auto& depo : depos) {
       std::cerr << " " << depo->id();
    }
    std::cerr << "\n";
}

void test_hydraTT()
{
    MyDepoSrc mds;
    MyHydraTT mh;
    const size_t ndepos = 10;
    for (size_t ind=0; ind<ndepos; ++ind) {

        MyHydraTT::input_queues iqs;
        size_t mode = ind % 4;
        if (mode % 2) {
            IDepo::pointer depo;
            mds(depo);
            std::get<1>(iqs).push_back(depo);
        }
        if (mode >= 2) {
            IDepo::pointer depo;
            mds(depo);
            std::get<0>(iqs).push_back(depo);
        }

        std::cerr << "input 0: "; dump_queue(std::get<0>(iqs));
        std::cerr << "input 1: "; dump_queue(std::get<1>(iqs));
        MyHydraTT::output_queues oqs;
        mh(iqs, oqs);
        std::cerr << "output: "; dump_queue(std::get<0>(oqs));
    }
}

void test_hydraVV()
{
    MyDepoSrc mds;
    MyHydraVV mh;
    const size_t ndepos = 10;
    for (size_t ind=0; ind<ndepos; ++ind) {

        MyHydraVV::input_queues iqs(2);
        size_t mode = ind % 4;
        if (mode % 2) {
            IDepo::pointer depo;
            mds(depo);
            iqs[1].push_back(depo);
        }
        if (mode >= 2) {
            IDepo::pointer depo;
            mds(depo);
            iqs[0].push_back(depo);
        }

        std::cerr << "input 0: "; dump_queue(iqs[0]);
        std::cerr << "input 1: "; dump_queue(iqs[1]);
        MyHydraVV::output_queues oqs(1);
        mh(iqs, oqs);
        std::cerr << "output: "; dump_queue(oqs[0]);
    }
}

void test_hydraVT()
{
    MyDepoSrc mds;
    MyHydraVT mh;
    const size_t ndepos = 10;
    for (size_t ind=0; ind<ndepos; ++ind) {

        MyHydraVT::input_queues iqs(2);
        size_t mode = ind % 4;
        if (mode % 2) {
            IDepo::pointer depo;
            mds(depo);
            iqs[1].push_back(depo);
        }
        if (mode >= 2) {
            IDepo::pointer depo;
            mds(depo);
            iqs[0].push_back(depo);
        }

        std::cerr << "input 0: "; dump_queue(iqs[0]);
        std::cerr << "input 1: "; dump_queue(iqs[1]);
        MyHydraVT::output_queues oqs;
        mh(iqs, oqs);
        std::cerr << "output: "; dump_queue(get<0>(oqs));
    }
}

void test_hydraTV()
{
    MyDepoSrc mds;
    MyHydraTV mh;
    const size_t ndepos = 10;
    for (size_t ind=0; ind<ndepos; ++ind) {

        MyHydraTV::input_queues iqs;
        size_t mode = ind % 4;
        if (mode % 2) {
            IDepo::pointer depo;
            mds(depo);
            std::get<1>(iqs).push_back(depo);
        }
        if (mode >= 2) {
            IDepo::pointer depo;
            mds(depo);
            std::get<0>(iqs).push_back(depo);
        }

        std::cerr << "input 0: "; dump_queue(std::get<0>(iqs));
        std::cerr << "input 1: "; dump_queue(std::get<1>(iqs));
        MyHydraTV::output_queues oqs(1);
        mh(iqs, oqs);
        std::cerr << "output: "; dump_queue(oqs[0]);
    }
}
int main()
{
    test_hydraTT();
    test_hydraVV();
    test_hydraVT();
    test_hydraTV();

    return 0;
}
