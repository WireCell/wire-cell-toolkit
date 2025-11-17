#include "WireCellSpng/CrossViewsExtract.h"
#include "WireCellSpng/SimpleTorchTensor.h"
#include "WireCellUtil/NamedFactory.h"

WIRECELL_FACTORY(SPNGCrossViewsExtract,
                 WireCell::SPNG::CrossViewsExtract,
                 WireCell::SPNG::CrossViewsExtract::fan_type,
                 WireCell::IConfigurable,
                 WireCell::INamed);

using WireCell::HanaJsonCPP::from_json;

using WireCell::HanaJsonCPP::to_json;


namespace WireCell::SPNG {

    CrossViewsExtract::CrossViewsExtract()
        : FanoutBase<ITorchTensor>("CrossViewsExtract")
    {
    }
    
    CrossViewsExtract::~CrossViewsExtract()
    {
    }

    void CrossViewsExtract::fanout_separate(const input_pointer& in, output_vector& outv)
    {
        logit(in, "cross views extract");
        Configuration md = in->metadata();
        auto ten = in->tensor();
        const size_t nout = outv.size();
        if (nout != m_ops.size()) {
            raise<ValueError>("mismatch between data (%d) and configured (%d) fanout size",
                              nout, m_ops.size());
        }

        for (size_t ind=0; ind<nout; ++ind) {
            auto op = m_ops[ind];
            auto out = op(ten);
            outv[ind] = std::make_shared<SimpleTorchTensor>(out, md);
        }
    }

    WireCell::Configuration CrossViewsExtract::default_configuration() const
    {
        // unusually, we go first to provide fan multiplicity.
        auto cfg = HanaJsonCPP::to_json(m_config);
        auto cfg2 = FanoutBase<ITorchTensor>::default_configuration();
        cfg2["multiplicity"] = (int)cfg["extraction"].size();
        update(cfg2, cfg);
        return cfg2;
    }


    void CrossViewsExtract::configure(const WireCell::Configuration& jconfig)
    {
        // unusually, we go first to provide fan multiplicity.
        from_json(m_config, jconfig);
        auto jcfg = jconfig;
        jcfg["multiplicity"] = (int)m_config.extraction.size();
        FanoutBase<ITorchTensor>::configure(jcfg);

        if (m_config.extraction.empty()) {
            raise<ValueError>("you probably did not mean to apply no extraction, check your config.");
        }

        // Build operations
        m_ops.clear();
        for (const auto& ex : m_config.extraction) {
            int mpn = 0;
            if (ex == "mp2") {
                mpn = 2;
            }
            else if (ex == "mp3") {
                mpn = 3;
            }
            else {
                raise<ValueError>("unsupported extraction: \"%s\"", ex);
            }
            const int mp = 1 << mpn;
            const int mask = (mp<<8)|(mp<<4)|mp;
            auto op = [mask](const torch::Tensor& tensor) {
                return torch::bitwise_and(tensor, mask).to(torch::kBool);
            };
            m_ops.push_back(op);
        }
    }


}


