/**
 * A simple implementation of IGeomService which provides a single
 */

#ifndef WIRECELLGEN_SIMPLEGEOMSERVICE
#define WIRECELLGEN_SIMPLEGEOMSERVICE

#include "WireCellIface/IGeomService.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellAux/Logger.h"
#include "WireCellUtil/Units.h"
#include "unordered_map"

namespace WireCell {
    namespace Clus {

        class SimpleGeomService : public Aux::Logger, public IGeomService, public IConfigurable {
           public:
            SimpleGeomService();
            virtual ~SimpleGeomService() {}

            // IConfigurable interface
            virtual WireCell::Configuration default_configuration() const;
            virtual void configure(const WireCell::Configuration& config);

            /// IGeomService interface
            virtual WireCell::Configuration get_params(const int apa, const int face) const;
            virtual bool is_in_FV(const WireCell::Point& point, const int apa, const int face) const;
            virtual bool is_in_FV_dim(const WireCell::Point& point, const int dim, const double margin, const int apa,
                                      const int face) const;
            virtual WireCell::Point get_corrected_point(const WireCell::Point& point,
                                                        const WireCell::IGeomService::CorrectionType type,
                                                        const int apa, const int face) const;

            struct FV {
                double m_FV_xmin{1 * units::cm};
                double m_FV_xmax{255 * units::cm};
                double m_FV_ymin{-99.5 * units::cm};
                double m_FV_ymax{101.5 * units::cm};
                double m_FV_zmin{15 * units::cm};
                double m_FV_zmax{1022 * units::cm};
                double m_FV_xmin_margin{2 * units::cm};
                double m_FV_xmax_margin{2 * units::cm};
                double m_FV_ymin_margin{2.5 * units::cm};
                double m_FV_ymax_margin{2.5 * units::cm};
                double m_FV_zmin_margin{3 * units::cm};
                double m_FV_zmax_margin{3 * units::cm};
            };
           private:
            WireCell::Configuration m_tpcparams;
            std::unordered_map<std::string, FV> m_FV_map;
        };
    }  // namespace Clus

}  // namespace WireCell

#endif
