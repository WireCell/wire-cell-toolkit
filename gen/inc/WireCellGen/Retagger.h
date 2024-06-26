/* The retagger support frame and trace tag rule sets like FrameFanin
 * and FrameFanout but also support a "merge" ruleset which can
 * combine trace sets.  For example an element of the tag_rules in
 * configuration may look like:
  {
    type: "Retagger",
    data: {
        tag_rules: [{
            merge: {
                "gauss\d": "gauss",
                "wiener\d": "wiener",
            }
        }],
     }
  }

 * This would put all trace sets tagged like gauss0, gauss1, etc into
 * an output trace set tagged "gauss".
 *
 * Note that although Retagger is a 1-1 frame filter, the tag_rules is
 * a list in order to match the data structure pattern with
 * FrameFanout and FrameFanin.
 *
 * This component, like others that do tag management, do not carry
 * forward tags by default. Any to be carried forward need to be given
 * in either the "frame" or "trace" rulesets.
 */

#ifndef WIRECELL_GEN_RETAGGER
#define WIRECELL_GEN_RETAGGER

#include "WireCellIface/IFrameFilter.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/TagRules.h"

namespace WireCell {

    namespace Gen {

        class Retagger : public IFrameFilter, public IConfigurable {
           public:
            Retagger();
            virtual ~Retagger();

            virtual void configure(const WireCell::Configuration& config);
            virtual WireCell::Configuration default_configuration() const;

            virtual bool operator()(const input_pointer& inframe, output_pointer& outframe);

           private:
            tagrules::Context m_trctx;
            std::map<std::string, std::string> m_maskmap;
        };
    }  // namespace Gen
}  // namespace WireCell

#endif
