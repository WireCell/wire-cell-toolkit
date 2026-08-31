/**
 * The idea is that there are objects of a base interface IComponent
 * which may also implement some other interface, here IPort.  The
 * objects are initially known through IComponent and a templated
 * function, transfer(), must be called on these objects but through a
 * middle sublcass which contains some needed type info, here Port<T>.
 *
 * The solution is to create a context (CallTransfer and
 * CallTransferT<> structs) which knows the type T and then register a
 * number of instances of this context using the type_info::name() of
 * the instantiated T.
 *
 * Then, the base interface IPort has a pure abstract method to return
 * the type_info::name() supplied by the subclass.  This is used to
 * look up the corresponding context which takes the base IPort
 * objects and dynamic_cast's them to the needed subclass.
 */

#include "WireCellUtil/doctest.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/Interface.h"
#include <typeinfo>
#include <map>
#include <string>

using namespace std;

struct IPort {
    virtual ~IPort(){};
    virtual std::string port_type_name() const = 0;
};

template <typename T>
struct Port : public IPort {
    typedef T port_type;
    virtual ~Port() {}
    virtual std::string port_type_name() const { return typeid(port_type).name(); }
    virtual bool put(const port_type& in) { return false; }
    virtual bool get(port_type& out) const { return false; }
    virtual T make() const { return 0; }
};

struct SubI : public Port<int> {
    virtual ~SubI() {}
    virtual bool get(port_type& out) const
    {
        out = 42;
        spdlog::debug("SubI::get({})", out);
        return true;
    }
};

struct SubF : public Port<float> {
    virtual ~SubF() {}
    virtual bool get(port_type& out) const
    {
        out = 6.9;
        spdlog::debug("SubF::get({})", out);
        return true;
    }
};

struct SubFIn : public Port<float> {
    virtual ~SubFIn() {}
    virtual bool put(const port_type& in)
    {
        spdlog::debug("SubFIn::put({})", in);
        return true;
    }
    virtual bool get(port_type& out) const
    {
        spdlog::debug("SubIn::get() can not provide data");
        return true;
    }
};

struct SubFOut : public Port<float> {
    virtual ~SubFOut() {}
    virtual bool put(const port_type& in)
    {
        spdlog::debug("SubFOut::put({}) can not accept data", in);
        return false;
    }
    virtual bool get(port_type& out) const
    {
        out = 6.9;
        spdlog::debug("SubFOut::get({})", out);
        return true;
    }
};

template <typename A, typename B>
bool transfer(const A& a, B& b)
{
    spdlog::debug("transfer: {}<->{} {}<->{}",
                  typeid(A).name(), typeid(a).name(),
                  typeid(B).name(), typeid(b).name());
    if (a.port_type_name() != b.port_type_name()) {
        spdlog::debug("Port type mismatch: {} != {}", a.port_type_name(), b.port_type_name());
        return false;
    }
    typename A::port_type dat;
    if (!a.get(dat)) {
        spdlog::debug("Failed to get output of type {}", a.port_type_name());
        return false;
    }
    if (!b.put(dat)) {
        spdlog::debug("Failed to put input of type {}", b.port_type_name());
        return false;
    }
    return true;
}

struct CallTransfer {
    virtual ~CallTransfer() {}
    virtual bool call(const IPort& a, IPort& b) = 0;
};

template <typename T>
struct CallTransferT : public CallTransfer {
    typedef T port_type;
    virtual ~CallTransferT() {}
    std::string port_type_name() const { return typeid(port_type).name(); }
    bool call(const IPort& a, IPort& b)
    {
        const Port<port_type>* pa = dynamic_cast<const Port<port_type>*>(&a);
        Port<port_type>* pb = dynamic_cast<Port<port_type>*>(&b);
        return transfer(*pa, *pb);
    }
};

static map<string, CallTransfer*> callers;

template <typename T>
void register_caller()
{
    CallTransferT<T>* ct = new CallTransferT<T>;
    callers[ct->port_type_name()] = ct;
}

CallTransfer* get_caller(const IPort& comp)
{
    const IPort* port = dynamic_cast<const IPort*>(&comp);
    if (!port) {
        return nullptr;
    }
    return callers[port->port_type_name()];
}

TEST_SUITE("typeid") {

TEST_CASE("typeid info") {
    IPort* si = new SubI;
    IPort* sf = new SubF;
    IPort* sfi = new SubFIn;
    IPort* sfo = new SubFOut;

    spdlog::debug("typeid(si) = {}", typeid(si).name());
    spdlog::debug("typeid(*si) = {}", typeid(*si).name());
    spdlog::debug("typeid(sf) = {}", typeid(sf).name());
    spdlog::debug("typeid(*sf) = {}", typeid(*sf).name());
    spdlog::debug("typeid(sfo) = {}", typeid(sfo).name());
    spdlog::debug("typeid(*sfo) = {}", typeid(*sfo).name());
    spdlog::debug("typeid(sfi) = {}", typeid(sfi).name());
    spdlog::debug("typeid(*sfi) = {}", typeid(*sfi).name());

    CHECK(si->port_type_name() == typeid(int).name());
    CHECK(sf->port_type_name() == typeid(float).name());
    CHECK(sfo->port_type_name() == typeid(float).name());
    CHECK(sfi->port_type_name() == typeid(float).name());

    delete si;
    delete sf;
    delete sfi;
    delete sfo;
}

TEST_CASE("dynamic cast transfer") {
    IPort* sfi = new SubFIn;
    IPort* sfo = new SubFOut;

    Port<float>* bfi = dynamic_cast<Port<float>*>(sfi);
    Port<float>* bfo = dynamic_cast<Port<float>*>(sfo);
    REQUIRE(bfi);
    REQUIRE(bfo);
    CHECK(transfer(*bfo, *bfi));

    delete sfi;
    delete sfo;
}

TEST_CASE("caller registry") {
    register_caller<int>();
    register_caller<float>();

    IPort* ci = new SubFIn;
    IPort* co = new SubFOut;

    CallTransfer* ct = get_caller(*ci);
    REQUIRE(ct);
    CHECK(ct->call(*co, *ci));

    delete ci;
    delete co;
}

}  // TEST_SUITE("typeid")
