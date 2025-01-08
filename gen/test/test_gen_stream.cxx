// test boost iostream with objects
// http://stackoverflow.com/questions/14086417/how-to-write-custom-input-stream-in-c

#include "WireCellUtil/Testing.h"

#include <string>
#include <memory>
#include <set>
#include <iostream>

struct MyData {
    int i;
    float f;
    double d;
    std::string s;
    MyData(int i = 0, float f = 0.0, double d = 0.0, const std::string& s = "")
      : i(i)
      , f(f)
      , d(d)
      , s(s)
    {
    }
};

std::ostream& operator<<(std::ostream& o, const MyData& md)
{
    o << "<MyData i:" << md.i << " f:" << md.f << " d:" << md.d << " s:\"" << md.s << "\">";
    return o;
}

typedef std::shared_ptr<const MyData> MyDataPtr;
typedef std::set<MyDataPtr> MyDataStore;

class IMySource {
   public:
    virtual ~IMySource() {}
    virtual MyDataPtr operator()() = 0;
    virtual bool eof() = 0;
};

class MyGen : virtual public IMySource {
    int m_count;
    int m_max;

   public:
    MyGen(int max = -1)
      : m_count(0)
      , m_max(max)
    {
    }

    virtual MyDataPtr operator()()
    {
        MyData* rawp = new MyData(m_count, m_count, m_count);
        ++m_count;
        return MyDataPtr(rawp);
    }
    virtual bool eof()
    {
        if (m_max < 0) return false;
        return m_count >= m_max;
    }
};

class MyFilter : virtual public IMySource {
    IMySource& m_src;
    int m_mod;

   public:
    MyFilter(IMySource& src, int m = 2)
      : m_src(src)
      , m_mod(m)
    {
    }
    virtual ~MyFilter(){};

    virtual MyDataPtr operator()()
    {
        while (!m_src.eof()) {
            MyDataPtr p = m_src();
            if (0 == p->i % m_mod) {
                return p;
            }
        }
        return 0;
    }
    virtual bool eof() { return m_src.eof(); }
};

using namespace std;

int main()
{
    MyGen mg(10);
    MyFilter mf(mg);
    MyDataPtr p;
    while ((p = mf())) {
        cerr << *p << endl;
    }
}
