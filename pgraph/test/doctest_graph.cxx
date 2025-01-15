#include "WireCellUtil/Logging.h"
#include "WireCellPgraph/Graph.h"
#include "WireCellUtil/doctest.h"
#include "WireCellUtil/String.h"

#include <boost/container_hash/hash.hpp>

#include <sstream>

using namespace WireCell;
using spdlog::debug;
using WireCell::String::format;

class IdNode : public Pgraph::Node {
   public:
    IdNode(const std::string& name, int id, size_t nin = 0, size_t nout = 0)
      : m_name(name)
      , m_id(id)
      , m_record{nullptr}
    {
        using Pgraph::Port;
        for (size_t ind = 0; ind < nin; ++ind) {
            m_ports[Port::input].push_back(Pgraph::Port(this, Pgraph::Port::input, "int"));
        }
        for (size_t ind = 0; ind < nout; ++ind) {
            m_ports[Port::output].push_back(Pgraph::Port(this, Pgraph::Port::output, "int"));
        }
    }
    int id() { return m_id; }

    virtual std::string ident()
    {
        std::stringstream ss;
        ss << m_name << "[" << m_id << "]";
        return ss.str();
    }
    void msg(const std::string s)
    {
        debug("#{} {}:\t{}", instance(), ident(), s);
    }

    virtual bool ready()
    {
        using Pgraph::Port;
        for (auto& p : m_ports[Port::input]) {
            if (p.empty()) return false;
        }
        return true;
    }

    void set_record(std::vector<int>* r) { m_record = r; }
    void record() { if (m_record) { m_record->push_back(m_id); } }

   private:
    std::string m_name;
    int m_id;
    std::vector<int>* m_record;
};
class Source : public IdNode {
   public:
    Source(int id, int beg, int end)
      : IdNode("src", id, 0, 1)
      , m_num(beg)
      , m_end(end)
    {
    }
    virtual ~Source() {}
    virtual bool ready() { return m_num < m_end; }
    virtual bool operator()()
    {
        record();
        if (m_num >= m_end) {
            msg("dry");
            return false;
        }
        msg(format("make: %d", m_num));
        Pgraph::Data d = m_num;
        oport().put(d);
        ++m_num;
        return true;
    }

   private:
    int m_num, m_end;
};
class Sink : public IdNode {
   public:
    Sink(int id)
      : IdNode("dst", id, 1, 0)
    {
    }
    virtual ~Sink() {}
    virtual bool operator()()
    {
        record();
        if (iport().empty()) {
            return false;
        }
        int d = boost::any_cast<int>(iport().get());
        msg(format("sink: %d", d));
        return true;
    }
};
class Njoin : public IdNode {
   public:
    Njoin(int id, int n)
      : IdNode("joi", id, n, 1)
    {
    }
    virtual bool ready()
    {
        auto& ip = input_ports();
        for (size_t ind = 0; ind < ip.size(); ++ind) {
            auto& p = ip[ind];
            if (p.empty()) {
                return false;
            }
        }
        return true;
    }
    virtual bool operator()()
    {
        record();
        Pgraph::Queue outv;
        std::stringstream ss;
        ss << "join: ";
        for (auto p : input_ports()) {
            if (p.empty()) {
                continue;
            }
            Pgraph::Data d = p.get();
            int n = boost::any_cast<int>(d);
            ss << n << " ";
            outv.push_back(d);
        }
        msg(ss.str());

        if (outv.empty()) {
            return false;
        }

        Pgraph::Data out = outv;
        oport().put(out);
        return true;
    }
};
class SplitQueueBuffer : public IdNode {
   public:
    SplitQueueBuffer(int id)
      : IdNode("sqb", id, 1, 1)
    {
    }
    virtual bool read()
    {
        if (!m_buf.empty()) {
            return true;
        }
        return IdNode::ready();
    }
    virtual bool operator()()
    {
        record();
        if (m_buf.empty()) {
            if (iport().empty()) {
                return false;
            }
            m_buf = boost::any_cast<Pgraph::Queue>(iport().get());
        }
        if (m_buf.empty()) {
            return false;
        }
        auto d = m_buf.front();
        m_buf.pop_front();
        oport().put(d);
        return true;
    }

   private:
    Pgraph::Queue m_buf;
};

class Nfan : public IdNode {
   public:
    Nfan(int id, int n)
      : IdNode("fan", id, 1, n)
    {
    }
    virtual bool operator()()
    {
        record();
        if (iport().empty()) {
            return false;
        }
        auto obj = iport().get();
        int d = boost::any_cast<int>(obj);
        msg(format("nfan: %d", d));
        for (auto p : output_ports()) {
            p.put(obj);
        }
        return true;
    }
};
class Func : public IdNode {
   public:
    Func(int id)
      : IdNode("fun", id, 1, 1)
    {
    }
    virtual bool operator()()
    {
        record();
        if (iport().empty()) {
            return false;
        }
        Pgraph::Data out = iport().get();
        int d = boost::any_cast<int>(out);
        msg(format("func: %d", d));
        oport().put(out);
        return true;
    }
};

// fixme: add N->M: hydra


using namespace WireCell;

static size_t do_graph(bool swap, bool extra)
{
    using Pgraph::Graph;
    using Pgraph::Node;

    int count = 0;
    // Hold shared pointers to mimic how WCT does it internally.
    // It's a bit silly to do it here.
    std::shared_ptr<IdNode> src1, src2;
    if (swap) {
        src2 = std::make_shared<Source>(count++, 0, 4);
        src1 = std::make_shared<Source>(count++, 10, 14);
    }
    else {
        src1 = std::make_shared<Source>(count++, 0, 4);
        src2 = std::make_shared<Source>(count++, 10, 14);
    }
    auto dst1 = std::make_shared<Sink>(count++);
    std::unique_ptr<int []> dummy;
    if (extra) {
        dummy = std::make_unique<int []>(1000);
    }
    auto dst2 = std::make_shared<Sink>(count++);
    auto fun1 = std::make_shared<Func>(count++);
    auto fun2 = std::make_shared<Func>(count++);
    auto fun3 = std::make_shared<Func>(count++);
    auto fan1 = std::make_shared<Nfan>(count++, 2);
    auto joi1 = std::make_shared<Njoin>(count++, 2);
    auto sqb1 = std::make_shared<SplitQueueBuffer>(count++);

    Graph g;
    g.connect(src1.get(), fun1.get());
    g.connect(fun1.get(), fan1.get());
    g.connect(fan1.get(), dst1.get());
    g.connect(fan1.get(), fun2.get(), 1);
    g.connect(fun2.get(), joi1.get());
    g.connect(src2.get(), fun3.get());
    g.connect(fun3.get(), joi1.get(), 0, 1);
    g.connect(joi1.get(), sqb1.get());
    g.connect(sqb1.get(), dst2.get());

    auto sorted = g.sort_kahn();
    debug("Sorted to {} nodes", count);

    std::vector<int> record;

    for (size_t ind = 0; ind < sorted.size(); ++ind) {
        IdNode* idn = dynamic_cast<IdNode*>(sorted[ind]);
        idn->msg(format("at index %d", ind));
        idn->set_record(&record);
    }
    debug("Executing:");

    g.execute();

    debug("{} executions", record.size());
    size_t checksum = 0;
    boost::hash_combine(checksum, record.size());
    for (int id : record) {
        boost::hash_combine(checksum, id);
    }
    debug("checksum: {}", checksum);
    return checksum;
}
TEST_CASE("pgraph original pipegraph test")
{
    auto cs1 = do_graph(false, false);
    auto cs2 = do_graph(true, false);
    auto cs3 = do_graph(false, true);
    auto cs4 = do_graph(true, true);

    // REQUIRE(checksum == 12973496632588361711ULL); // plant a potential bomb
    // 16328357942837688713 == 12973496632588361711

    REQUIRE(cs1 != cs2);
    REQUIRE(cs3 != cs4);
    REQUIRE(cs1 == cs3);
    REQUIRE(cs2 == cs4);
}
