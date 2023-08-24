#include "WireCellSio/TensorFileSource.h"
#include "WireCellUtil/Stream.h"
#include "WireCellUtil/String.h"
#include "WireCellUtil/Dtype.h"
#include "WireCellAux/SimpleTensorSet.h"

#include "WireCellUtil/NamedFactory.h"

WIRECELL_FACTORY(TensorFileSource, WireCell::Sio::TensorFileSource,
                 WireCell::INamed,
                 WireCell::ITensorSetSource,
                 WireCell::IConfigurable)


using namespace WireCell;
using namespace WireCell::Sio;
using namespace WireCell::Aux;
using namespace WireCell::Stream;
using namespace WireCell::String;


TensorFileSource::TensorFileSource()
    : Aux::Logger("TensorFileSource", "sio")
{
}

TensorFileSource::~TensorFileSource()
{
}


WireCell::Configuration TensorFileSource::default_configuration() const
{
    Configuration cfg;
    cfg["inname"] = m_inname;
    cfg["prefix"] = m_prefix;
    return cfg;
}

void TensorFileSource::configure(const WireCell::Configuration& cfg)
{
    m_inname = get(cfg, "inname", m_inname);
    m_prefix = get(cfg, "prefix", m_prefix);

    m_in.clear();
    input_filters(m_in, m_inname);
    if (m_in.empty()) {
        THROW(ValueError() << errmsg{"TensorFileSource: unsupported inname: " + m_inname});
    }

    log->debug("reading file={} with prefix={}", m_inname, m_prefix);
}

void TensorFileSource::finalize()
{
}

/*
  <prefix>tensorset_<ident>_metadata.json 
  <prefix>tensor_<ident>_<index>_metadata.npy
  <prefix>tensor_<ident>_<index>_array.json
*/
struct ParsedFilename {
    enum Type { bad, set, ten };
    enum Form { unknown, npy, json };
    Type type{bad};
    Form form{unknown};
    int ident{-1};
    int index{-1};              // only for ten
};
static
ParsedFilename parse_fname(const std::string& fname,
                           const std::string& prefix)
{
    ParsedFilename ret;
    if (! startswith(fname, prefix)) {
        return ret;
    }

    std::string basename;
    if (endswith(fname, ".json")) {
        ret.form = ParsedFilename::json;
        basename = fname.substr(0, fname.size() - 5);
    }
    else if (endswith(fname, ".npy")) {
        ret.form = ParsedFilename::npy;
        basename = fname.substr(0, fname.size() - 4);
    }
    else return ret;

    auto parts = split(basename, "_");
    if (parts.size() == 3 and endswith(parts[0], "tensorset") and parts[2] == "metadata") {
        ret.type = ParsedFilename::set;
        ret.ident = atoi(parts[1].c_str());
        return ret;
    }
    else if (parts.size() == 4 and endswith(parts[0], "tensor")) {
        if (parts[3] == "metadata") {
            ret.type = ParsedFilename::ten;
            ret.form = ParsedFilename::json;
            ret.ident = atoi(parts[1].c_str());
            ret.index = atoi(parts[2].c_str());
            return ret;
        }
        else if (parts[3] == "array") {
            ret.type = ParsedFilename::ten;
            ret.form = ParsedFilename::npy;
            ret.ident = atoi(parts[1].c_str());
            ret.index = atoi(parts[2].c_str());
            return ret;
        }
        else return ret;
    }
    else return ret;
}

static
Configuration load_json(std::istream& in, size_t fsize)
{
    Configuration cfg;

    // Read into a buffer string to assure we consume fsize bytes and
    // protect against malformed JSON or leading/trailing spaces.
    std::string buffer;
    buffer.resize(fsize);
    in.read(buffer.data(), buffer.size());
    if (!in) { return cfg; }

    std::istringstream ss(buffer);
    ss >> cfg;
    return cfg;
}

struct TFSTensor : public WireCell::ITensor {
    // pigenc::File pig;
    TFSTensor() {}
    virtual ~TFSTensor() {}
    // explicit TFSTensor(pigenc::File&& pig)
    //     : pig(std::move(pig))
    // {
    // }
    size_t m_element_size{0};
    std::vector<std::byte> m_store;
    std::string m_dtype{""};
    shape_t m_shape;
    Configuration m_cfg;
    explicit TFSTensor(const pigenc::File& pig, Configuration cfg=Json::objectValue)
        : m_element_size(pig.header().type_size())
        , m_dtype(pig.header().dtype())
        , m_shape(pig.header().shape())
        , m_cfg(cfg)
    {
        const auto& pdata = pig.data();
        const size_t size = pdata.size();
        const std::byte* data = reinterpret_cast<const std::byte*>(pdata.data());
        m_store.insert(m_store.end(), data, data + size);
    }

    // An ITensor may not have an array part.
    explicit TFSTensor(Configuration cfg)
        : m_cfg(cfg)
    {
    }
    virtual const std::type_info& element_type() const
    {
        // return dtype_info(dtype());
        return dtype_info(m_dtype);
    }
    virtual size_t element_size() const
    {
        // return pig.header().type_size();
        return m_element_size;
    }
    virtual std::string dtype() const
    {
        // return pig.header().dtype();
        return m_dtype;
    }
    virtual shape_t shape() const
    {
        // return pig.header().shape();
        return m_shape;
    }
    virtual const std::byte* data() const
    {
        // return (const std::byte*)pig.data().data();
        return m_store.data();
    }
    virtual size_t size() const
    {
        // return pig.header().data_size();
        return m_store.size();
    }
    virtual Configuration metadata() const
    {
        return m_cfg; 
    }
};

ITensorSet::pointer TensorFileSource::load()
{
    int ident = -1;

    Configuration setmd;

    struct TenInfo {
        Configuration md;
        std::shared_ptr<TFSTensor> ten;
    };
    std::map<int, TenInfo> teninfo; // by index

    while (true) {

        if (m_cur.fsize == 0) {
            clear();
            custard::read(m_in, m_cur.fname, m_cur.fsize);
            if (m_in.eof()) {
                // log->debug("call={}, read stream EOF from file={}",
                //            m_count, m_inname);
                break;
            }
            if (!m_in) {
                log->critical("call={}, read stream error with file={}",
                              m_count, m_inname);
                return nullptr;
            }
            if (!m_cur.fsize) {
                log->critical("call={}, short read from file={}",
                              m_count, m_inname);
                return nullptr;
            }
        }

        auto pf = parse_fname(m_cur.fname, m_prefix);

        // log->debug("read file={} size={} type={} form={} ident={} index={}",
        //            m_cur.fname, m_cur.fsize,
        //            pf.type, pf.form, pf.ident, pf.index);

        if (pf.type == ParsedFilename::bad or pf.form == ParsedFilename::unknown) {
            m_in.seekg(m_cur.fsize, m_in.cur);
            clear();
            continue;
        }
        if (ident < 0) {
            ident = pf.ident;   // first time through
        }
        if (ident != pf.ident) {
            break;              // started reading next tensor set.
        }

        // tensor set md 
        if (pf.type == ParsedFilename::set) {
            setmd = load_json(m_in, m_cur.fsize);
            clear();
            continue;
        }

        // tensor md or array
        if (pf.type == ParsedFilename::ten) {
            if (pf.form == ParsedFilename::json) {
                teninfo[pf.index].md = load_json(m_in, m_cur.fsize);
                clear();
            }
            else {
                pigenc::File pig;
                pig.read(m_in);
                teninfo[pf.index].ten = std::make_shared<TFSTensor>(std::move(pig));
                clear();
            }
            continue;
        }
    }
    if (ident < 0) {
        return nullptr;
    }

    auto sv = std::make_shared<ITensor::vector>();
    for (auto& [ind, ti] : teninfo) {
        if (!ti.ten) {          // an array-less tensor
            ti.ten = std::make_shared<TFSTensor>(ti.md);
        }
        else {
            ti.ten->m_cfg = ti.md;
        }
        sv->push_back(ti.ten);
    }
    return std::make_shared<SimpleTensorSet>(ident, setmd, sv);
}

void TensorFileSource::clear()
{
    m_cur = header_t();
}

bool TensorFileSource::operator()(ITensorSet::pointer &out)
{
    out = nullptr;
    if (m_eos_sent) {
        log->debug("past EOS at call={}", m_count++);
        return false;
    }
    out = load();
    if (!out) {
        m_eos_sent = true;
    }
    if (!out) {
        log->debug("sending EOS at call={}", m_count++);
    }
    else {
        log->debug("read done at call={}", m_count++);
    }
    return true;
    
}

