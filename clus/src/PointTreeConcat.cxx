// PointTreeConcat -- join N tensor sets that carry point-cloud trees at
// DIFFERENT datapaths into one tensor set (doc pdvd/125).
//
// PointTreeMerging merges the "live" (and "dead") trees of its inputs into ONE
// tree.  This node does not merge: it concatenates the tensors as they are, so
// a job can hand MultiAlgBlobClustering a second, separate grouping (for
// example the deghosted cells next to the production "live" tree; MABC then
// lists it in "groupings" and maps its datapath with "insubpaths").  The ident
// and the metadata of the output are those of input 0.  The caller is
// responsible for the datapaths being distinct.  A NEW component: absent from
// every existing graph => no existing output changes.

#include "WireCellAux/Logger.h"
#include "WireCellAux/SimpleTensorSet.h"

#include "WireCellIface/ITensorSetFanin.h"
#include "WireCellIface/IConfigurable.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Exceptions.h"

namespace WireCell::Clus {

    class PointTreeConcat : public Aux::Logger, public ITensorSetFanin, public IConfigurable {
       public:
        PointTreeConcat()
          : Aux::Logger("PointTreeConcat", "clus")
        {
        }
        virtual ~PointTreeConcat() = default;

        virtual void configure(const WireCell::Configuration& cfg)
        {
            m_multiplicity = get<int>(cfg, "multiplicity", m_multiplicity);
            if (m_multiplicity < 1) {
                raise<ValueError>("PointTreeConcat: multiplicity must be positive");
            }
        }
        virtual WireCell::Configuration default_configuration() const
        {
            Configuration cfg;
            cfg["multiplicity"] = m_multiplicity;
            return cfg;
        }

        // INode, override because we get multiplicity at run time.
        virtual std::vector<std::string> input_types()
        {
            const std::string tname = std::string(typeid(input_type).name());
            return std::vector<std::string>(m_multiplicity, tname);
        }

        virtual bool operator()(const input_vector& invec, output_pointer& out)
        {
            out = nullptr;
            if ((int) invec.size() != m_multiplicity) {
                raise<ValueError>("PointTreeConcat: unexpected multiplicity got %d want %d", invec.size(), m_multiplicity);
            }
            size_t neos = 0;
            for (const auto& in : invec) {
                if (!in) { ++neos; }
            }
            if (neos == invec.size()) {
                SPDLOG_LOGGER_DEBUG(log, "EOS at call {}", m_count++);
                return true;
            }
            if (neos) {
                raise<ValueError>("PointTreeConcat: %d of %d inputs at EOS", neos, invec.size());
            }
            auto tens = std::make_shared<ITensor::vector>();
            for (size_t i = 0; i < invec.size(); ++i) {
                const auto its = invec[i]->tensors();
                SPDLOG_LOGGER_DEBUG(log, "input[{}] ident={} tensors={}", i, invec[i]->ident(), its ? its->size() : 0);
                if (its) tens->insert(tens->end(), its->begin(), its->end());
            }
            out = std::make_shared<Aux::SimpleTensorSet>(invec[0]->ident(), invec[0]->metadata(), tens);
            ++m_count;
            return true;
        }

       private:
        int m_multiplicity{2};
        size_t m_count{0};
    };

}  // namespace WireCell::Clus

WIRECELL_FACTORY(PointTreeConcat, WireCell::Clus::PointTreeConcat,
                 WireCell::INamed,
                 WireCell::ITensorSetFanin,
                 WireCell::IConfigurable)
