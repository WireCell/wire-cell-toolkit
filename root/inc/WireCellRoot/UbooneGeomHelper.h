/**
 * UBoone-specific implementation of IClusGeomHelper that provides
 * data-driven Space Charge Effect (SCE) position correction.
 *
 * The correction is loaded from a ROOT TH3F file
 * (SCEoffsets_dataDriven_combined_bkwd_Jan18.root) containing histograms
 * "hDx", "hDy", "hDz".  Logic ported from WCP::TPCParams::func_pos_SCE_correction
 * in prototype_base/data/src/TPCParams.cxx.
 *
 * Lives in WireCellRoot (not WireCellClus) because it depends on ROOT.
 */

#ifndef WIRECELLROOT_UBOONEGEOMHELPER
#define WIRECELLROOT_UBOONEGEOMHELPER

#include "WireCellClus/IClusGeomHelper.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellAux/Logger.h"
#include "WireCellUtil/Units.h"
#include <unordered_map>
#include <string>

// ROOT forward declarations
class TH3D;

namespace WireCell {
    namespace Root {

        class UbooneGeomHelper : public Aux::Logger,
                                 public IClusGeomHelper,
                                 public IConfigurable {
           public:
            UbooneGeomHelper();
            virtual ~UbooneGeomHelper();

            // IConfigurable interface
            virtual WireCell::Configuration default_configuration() const;
            virtual void configure(const WireCell::Configuration& config);

            // IClusGeomHelper interface
            virtual WireCell::Configuration get_params(const int apa, const int face) const;
            virtual bool is_in_FV(const WireCell::Point& point, const int apa, const int face) const;
            virtual bool is_in_FV_dim(const WireCell::Point& point, const int dim, const double margin,
                                      const int apa, const int face) const;
            virtual WireCell::Point get_corrected_point(const WireCell::Point& point,
                                                        const WireCell::IClusGeomHelper::CorrectionType type,
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

            // SCE position correction histograms (UBoone data-driven).
            // Non-null only after a valid sce_offsets_file is configured.
            TH3D* m_h3_Dx{nullptr};
            TH3D* m_h3_Dy{nullptr};
            TH3D* m_h3_Dz{nullptr};

            void load_sce_offsets(const std::string& filename);
        };

    }  // namespace Root
}  // namespace WireCell

#endif
