#include "WireCellSpng/Testing.h"
#include "WireCellSpng/HanaConfigurable.h"
#include "WireCellUtil/Configuration.h"
#include "WireCellIface/IConfigurable.h"

using namespace WireCell;


// "simple" configurable
struct SC : public virtual IConfigurable {
    virtual void configure(const WireCell::Configuration& jconfig) {
        std::cerr << "SC::configure: " << jconfig << "\n";
    }
    virtual WireCell::Configuration default_configuration() const {
        Configuration cfg;
        cfg["name"] = "SC";
        cfg["sc_name"] = "SC";
        std::cerr << "SC::default_configuration: " << cfg << "\n";
        return cfg;
    }
};
struct SC2 : public virtual IConfigurable {
    virtual void configure(const WireCell::Configuration& jconfig) {
        std::cerr << "SC2::configure: " << jconfig << "\n";
    }
    virtual WireCell::Configuration default_configuration() const {
        Configuration cfg;
        cfg["name"] = "SC2";
        cfg["sc2_name"] = "SC2";
        std::cerr << "SC2::default_configuration: " << cfg << "\n";
        return cfg;
    }
};

struct SCSC : public SC, public SC2,  public virtual IConfigurable {
    virtual void configure(const WireCell::Configuration& jconfig) {
        this->SC::configure(jconfig);
        this->SC2::configure(jconfig);
        std::cerr << "SCSC::configure: " << jconfig << "\n";
    }
    virtual WireCell::Configuration default_configuration() const {
        auto cfg = this->SC2::default_configuration();
        auto cfg2 = this->SC::default_configuration();
        update(cfg, cfg2);
        cfg["name"] = "SCSC";
        cfg["scsc_name"] = "SCSC";
        std::cerr << "SCSC::default_configuration: " << cfg << "\n";
        return cfg;
    }
};

// Here is example of aggregation by only member variable.
// Actually, you can inherit from SC here but not include "this" in the dispatch.
struct SCwithSC /*: public SC, public virtual IConfigurable*/ {
    SC sc;
    SC2 sc2;

    WireCell::Configuration default_configuration() const {
        return default_configuration_each({&sc, &sc2});
    }
    virtual void configure(const WireCell::Configuration& jconfig) {
        configure_each({&sc, &sc2}, jconfig);
    }
    
};

struct AC1 : AutoConfigurable<AC1>
{
    virtual void configured() {
        std::cerr << "AC1::configure_json: " << m_json_config << "\n";
    }
    virtual WireCell::Configuration default_configuration_json() const {
        Configuration cfg;
        cfg["name"] = "AC1";
        cfg["ac_name"] = "AC1";
        std::cerr << "AC1::default_configuration_json: " << cfg << "\n";
        return cfg;
    }
};

struct ACSC : AutoConfigurable<ACSC, SC>
{
    
    virtual void configured() {
        std::cerr << "ACSC::configure_json: " << m_json_config << "\n";
    }

    virtual WireCell::Configuration default_configuration_json() const {
        Configuration cfg;
        cfg["name"] = "ACSC";
        cfg["ac_name"] = "ACSC";
        std::cerr << "ACSC::default_configuration_json: " << cfg << "\n";
        return cfg;

    }
};

// Three configs with unique and potentially conflicting attributes.
struct HCC1 {
    int number = 41;
    std::string moto_hcc1 = "I'm HCC1";
};
BOOST_HANA_ADAPT_STRUCT(HCC1, number, moto_hcc1);
struct HCC2 {
    int number = 42;
    std::string moto_hcc2 = "I'm HCC2";
};
BOOST_HANA_ADAPT_STRUCT(HCC2, number, moto_hcc2);
struct HCC3 {
    int number = 43;
    std::string moto_hcc3 = "I'm HCC3";
};
BOOST_HANA_ADAPT_STRUCT(HCC3, number, moto_hcc3);

struct HCSC1 : HanaConfigurable<HCSC1, HCC1, SC> {
    virtual void configured() {
        std::cerr << "HCSC1 configured: number=" << m_hana_config.number << " moto=\"" << m_hana_config.moto_hcc1 << "\"\n";
    }
};

struct HCAC1 : HanaConfigurable<HCAC1, HCC1, ACSC> {
    virtual void configured() {
        std::cerr << "HCAC1 configured: number=" << m_hana_config.number << " moto=\"" << m_hana_config.moto_hcc1 << "\"\n";
    }
};
struct HCSC2 : HanaConfigurable<HCSC2, HCC2, SC> {
    virtual void configured() {
        std::cerr << "HCSC2 configured: number=" << m_hana_config.number << " moto=\"" << m_hana_config.moto_hcc2 << "\"\n";
    }
};

struct HCAC2 : HanaConfigurable<HCAC2, HCC2, ACSC> {
    virtual void configured() {
        std::cerr << "HCAC2 configured: number=" << m_hana_config.number << " moto=\"" << m_hana_config.moto_hcc2 << "\"\n";
    }
};

struct HCHC : HanaConfigurable<HCHC, HCC3, HCSC1, HCSC2> {
    virtual void configured() {
        std::cerr << "HCHC configured: number=" << m_hana_config.number << " moto=\"" << m_hana_config.moto_hcc3 << "\"\n";
    }
};

struct HConly : HanaConfigurable<HConly, HCC3> {
    virtual void configured() {
        std::cerr << "HConly configured: number=" << m_hana_config.number << " moto=\"" << m_hana_config.moto_hcc3 << "\"\n";
    }
};

struct HCwithHC : public HanaConfigurable<HCwithHC, HCC1, HConly>, public virtual IConfigurable {
    HCAC2 hcac;
    HCSC2 hcsc;


    WireCell::Configuration default_configuration_extra() const {
        return default_configuration_each({&hcac, &hcsc});
    }
    virtual void configure_extra(const WireCell::Configuration& jconfig) {
        configure_each({&hcac, &hcsc}, jconfig);
    }

    virtual void configured() {
        std::cerr << "HCwithHC configured: number=" << m_hana_config.number << " moto=\"" << m_hana_config.moto_hcc1 << "\"\n";
    }

    
};

// Doesn't work because SC is also base of ACSC
// warning: direct base ‘SC’ inaccessible in ‘WireCell::HanaConfigurable<HCACSC, HCC, AC, SC>’ due to ambiguity [-Winaccessible-base]
// struct HCACSC : HanaConfigurable<HCACSC, HCC, ACSC, SC> {
//     virtual void configured() {
//         std::cerr << "HCACSC configured: number=" << m_hana_config.number << " moto=\"" << m_hana_config.moto << "\"\n";
//     }
// };


template<typename C>
void cycle(C& c)
{
    std::cerr << "cycle: getting default config:\n";
    auto cfg = c.default_configuration();
    std::cerr << "cycle: default: config" << cfg << "\ncycle: configuring:\n";
    c.configure(cfg);
    std::cerr << "cycle: done\n";
}

TEST_CASE("spng simple configurable")
{
    SC sc;
    cycle(sc);
}
TEST_CASE("spng simple configurable inheritance")
{
    SCSC sc;
    cycle(sc);
}
TEST_CASE("spng auto configurable")
{
    ACSC ac;
    cycle(ac);
}
TEST_CASE("spng auto+simple configurable")
{
    ACSC ac;
    cycle(ac);
}
TEST_CASE("spng hana configurable single")
{
    HCSC1 hcsc;
    cycle(hcsc);
    HCAC1 hcac;
    cycle(hcac);
    // HCACSC hcacsc;
    // cycle(hcacsc);
}

TEST_CASE("spng hana configurable double")
{
    // this does not configure the HCSC1 and HCSC2 parents.
    HCHC h;
    cycle(h);
}

struct HC : HanaConfigurable<HC, HCC1> {
    virtual void configured() {
        std::cerr << "HC configured: number=" << m_hana_config.number << " moto=\"" << m_hana_config.moto_hcc1 << "\"\n";
    }
};
struct HC2 : HanaConfigurable<HC2, HCC2, HC> {
    virtual void configured() {
        std::cerr << "HC2 configured: number=" << m_hana_config.number << " moto=\"" << m_hana_config.moto_hcc2 << "\"\n";
    }
};

TEST_CASE("spng hana configurable simple and inherited")
{
    HC hc;
    cycle(hc);
    HC2 hc2;
    cycle(hc2);

}

TEST_CASE("spng simple configurable with member")
{
    SCwithSC c;
    cycle(c);
}
TEST_CASE("spng hana configurable with member")
{
    HCwithHC c;
    cycle(c);
}

TEST_CASE("spng hana helper functions")
{
    SC sc;
    SC2 sc2;
    AC1 ac;
    HConly hc;
    auto cfg = default_configuration_types(sc, sc2, ac, hc);
    std::cout << "cfg from many types: " << cfg << "\n";
    configure_types(cfg, sc, sc2, ac, hc);
        
}
