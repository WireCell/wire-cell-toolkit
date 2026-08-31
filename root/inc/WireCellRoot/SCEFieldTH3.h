/** ISCEField implementation backed by per-TPC TH3F histograms in a
 *  ROOT file (SBND dualmap layout: TrueBkwd_Displacement_{X,Y,Z}_{E,W}).
 *  IConfigurable so jsonnet wiring works the usual way; the actual
 *  SCE-correction transform in clus/ receives an ISCEField::pointer and
 *  never sees ROOT.
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
            virtual double displacement_y(int apa, double x, double y, double z) const override;
            virtual double displacement_z(int apa, double x, double y, double z) const override;

            // IConfigurable
            virtual void configure(const WireCell::Configuration& cfg) override;
            virtual WireCell::Configuration default_configuration() const override;

           private:
            // configuration
            std::string m_sce_map_file;
            std::string m_th3_name_E;     // X (mandatory)
            std::string m_th3_name_W;
            std::string m_th3_name_E_y;   // Y (optional transverse)
            std::string m_th3_name_W_y;
            std::string m_th3_name_E_z;   // Z (optional transverse)
            std::string m_th3_name_W_z;
            bool   m_axis_unit_mm{false};
            double m_sign{-1.0};

            // loaded state (all owned by m_file)
            std::unique_ptr<TFile> m_file;
            TH3F* m_hE{nullptr};
            TH3F* m_hW{nullptr};
            TH3F* m_hE_y{nullptr};
            TH3F* m_hW_y{nullptr};
            TH3F* m_hE_z{nullptr};
            TH3F* m_hW_z{nullptr};

            Log::logptr_t l;

            // Shared: pick E/W histo by apa, clamp, convert units, apply sign.
            // Returns 0 if the requested histo is absent (no-op component).
            double interp(TH3F* hE, TH3F* hW, int apa,
                          double x, double y, double z) const;

            static double clamp_axis(const TAxis* a, double v);

            // Load one (E,W) histo pair; hard=true throws on missing,
            // hard=false warns and leaves pointers null.
            void load_pair(const std::string& nameE, const std::string& nameW,
                           TH3F*& hE, TH3F*& hW, bool hard);
        };

    }  // namespace Root
}  // namespace WireCell

#endif
