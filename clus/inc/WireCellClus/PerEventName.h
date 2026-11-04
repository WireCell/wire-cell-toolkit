/** Per-event output names for a process that runs a group of events.
 *
 * The writers that make one file per event (the tracking visitors, the flat
 * tree, the tagger output, PrDisplayDump, the MABC's own Bee zip) take a name
 * TEMPLATE (clus.jsonnet evt_subdir, wct-pr.jsonnet).  Historically it was a
 * boost::format of the set ident alone ('pr_evt%1%', the art event number),
 * which collides when one process holds several subruns whose event numbers
 * restart (the gen2 BNB production: 20 subruns x events 1..10 per lar job).
 * The tokens below add the resolved run/subrun/event.
 */
#ifndef WIRECELLCLUS_PEREVENTNAME
#define WIRECELLCLUS_PEREVENTNAME

#include <string>

namespace WireCell::Clus {
    namespace Facade { class Ensemble; }

    /// Expand a per-event name template: the tokens {run}, {subrun}, {event}
    /// and {ident} are substituted, then, if a '%' remains, boost::format
    /// with the ident ('pr_evt%1%', the original idiom).  A template with
    /// neither is returned unchanged, so a literal name stays literal.
    std::string per_event_name(const std::string& tmpl, int ident, int run, int subrun, int event);

    /// The same from an ensemble: ident = ensemble.ident(); the RSE is the
    /// triplet the group's MultiAlgBlobClustering resolved (rse_valid()),
    /// else the fallback, -1 by default.  Writers that share one per-event file
    /// (the tracking visitor RECREATEs it, the flat tree UPDATEs it) must use
    /// the SAME fallback, so every caller uses the default (PR 544 review, item 10).
    std::string per_event_name(const std::string& tmpl, const Facade::Ensemble& ensemble,
                               int run = -1, int subrun = -1, int event = -1);

    /// True when the template names a different file per event: a '%', or one
    /// of the tokens above (any other '{' is literal).
    bool per_event_templated(const std::string& tmpl);
}

#endif
