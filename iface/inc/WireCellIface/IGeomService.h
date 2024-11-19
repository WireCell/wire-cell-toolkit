/**
 * Interface for a service which provides geometry information.
 */

#ifndef WIRECELLIFACES_IGEOMSERVICE
#define WIRECELLIFACES_IGEOMSERVICE

#include "WireCellUtil/IComponent.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/Point.h"

namespace WireCell {

    class IGeomService : public IComponent<IGeomService> {
       public:
        virtual ~IGeomService();

        /// Return WCP style TPC parameters in json format
        virtual Configuration get_params(const int apa, const int face) const = 0;

        virtual bool is_in_FV (const Point& point, const int apa, const int face) const = 0;

        /// true if pos_{dim} >= min_{dim} + margin and pos_{dim} <= max_{dim} - margin
        virtual bool is_in_FV_dim (const Point& point, const int dim, const double margin, const int apa, const int face) const = 0;

        enum CorrectionType {NONE=-1, SCE};
        virtual Point get_corrected_point(const Point& point, const CorrectionType type, const int apa, const int face) const = 0;

    };

}  // namespace WireCell

#endif
