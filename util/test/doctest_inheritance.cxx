#include "WireCellUtil/doctest.h"

#include <string>
#include <iostream>
#include <memory>

struct IBase {
    virtual ~IBase() = default;
};
struct IName : virtual public IBase {
    virtual std::string name() = 0;
    virtual ~IName() = default;
};
struct A : virtual public IName {
    virtual std::string name() { return "A"; }
    virtual ~A() = default;
};
struct B : public A {
    virtual std::string name() { return "B" + this->A::name(); }
    virtual ~B() = default;
};
struct C : public A {
    virtual std::string name() { return "C" + this->A::name(); }
    virtual ~C() = default;
};
struct D : public C, public B {
    virtual std::string name() { return "D" + this->C::name() + this->B::name(); }
    virtual ~D() = default;
};

TEST_CASE("diamond inheritance") {

  std::shared_ptr<IBase> ibase = std::make_shared<D>();
  std::shared_ptr<IName> id = std::dynamic_pointer_cast<IName>(ibase);
  std::cerr << id->name() << "\n";
  CHECK("DCABA" == id->name());
}
