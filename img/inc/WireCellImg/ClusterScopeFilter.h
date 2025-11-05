/**
 * ClusterScopeFilter
 */
#ifndef WIRECELL_CLUSTERSCOPEFILTER_H
#define WIRECELL_CLUSTERSCOPEFILTER_H

#include "WireCellIface/IClusterFilter.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellAux/Logger.h"

namespace WireCell {

    namespace Img {

        class ClusterScopeFilter : public Aux::Logger, public IClusterFilter, public IConfigurable {
           public:
            /// TODO: what is needed here
            /// FIXME: bit operation would be better
            enum BLOB_QUALITY_BITPOS { GOOD, BAD, POTENTIAL_GOOD, POTENTIAL_BAD, TO_BE_REMOVED };

            using vertex_tags_t = std::unordered_map<cluster_vertex_t, int>;
            // using vertex_tagging_t = std::function<void(const cluster_graph_t&, vertex_tags_t&)>;

            ClusterScopeFilter();
            virtual ~ClusterScopeFilter();

            virtual void configure(const WireCell::Configuration& cfg);
            virtual WireCell::Configuration default_configuration() const;

            virtual bool operator()(const input_pointer& in, output_pointer& out);

           private:
            int m_face_index{-1};
        };

    }  // namespace Img

}  // namespace WireCell

#endif /* WIRECELL_CLUSTERSCOPEFILTER_H */
