/** ISCEField implementation backed by per-TPC TH3F histograms in a
 *  ROOT file (SBND dualmap layout: TrueBkwd_Displacement_X_E,
 *  TrueBkwd_Displacement_X_W).  IConfigurable so jsonnet wiring works
 *  the usual way; the actual SCE-correction transform in clus/
 *  receives an ISCEField::pointer and never sees ROOT.
 */

#ifndef WIRECELLROOT_SCEFIELDTH3
#define WIRECELLROOT_SCEFIELDTH3

#include "WireCellIface/ISCEField.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/Logging.h"

#include <memory>
#include <string>

class TFile;
class TH3F;
class TAxis;

namespace WireCell {
    namespace Root {

        class SCEFieldTH3 : public WireCell::ISCEField,
                            public WireCell::IConfigurable {
           public:
            SCEFieldTH3();
            virtual ~SCEFieldTH3();

            // ISCEField
            virtual double displacement_x(int apa, double x, double y, double z) const override;

            // IConfigurable
            virtual void configure(const WireCell::Configuration& cfg) override;
            virtual WireCell::Configuration default_configuration() const override;

           private:
            // configuration
            std::string m_sce_map_file;
            std::string m_th3_name_E;
            std::string m_th3_name_W;
            bool   m_axis_unit_mm{false};
            double m_sign{-1.0};

            // loaded state
            std::unique_ptr<TFile> m_file;
            TH3F* m_hE{nullptr};  // owned by m_file
            TH3F* m_hW{nullptr};

            Log::logptr_t l;

            static double clamp_axis(const TAxis* a, double v);
        };

    }  // namespace Root
}  // namespace WireCell

#endif
