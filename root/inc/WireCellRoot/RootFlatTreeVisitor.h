/** RootFlatTreeVisitor: every tree of a per-event tracking ROOT file as ONE
 * per-event entry of a flat, CAF-style tree (ai-helper issue 39).
 *
 * Runs as the LAST IEnsembleVisitor of the PR pipeline, after every writer of
 * the per-event file (SbndPrMagnifyTrackingVisitor, UbooneTaggerOutputVisitor).
 * It reopens that file (UPDATE), reads every TTree in it and appends one tree,
 * `tree_name` (default "recTreeWireCell"), with exactly ONE entry for the event.
 * The flattening follows the flat CAF / SRProxy convention, so the event's
 * content can be appended entry by entry to a flat CAF's recTree:
 *
 *   <p>.run, <p>.subrun, <p>.evt   Int_t, from Trun (runNo/subRunNo/eventNo)
 *   <p>.valid                      Int_t, 1 (a merger writes 0 for an event
 *                                  without Wire-Cell output)
 *   <p>.<T>..length                Int_t, the number of rows of tree T
 *   a scalar branch b of T         <p>.<T>.<b>[<p>.<T>..length], same type
 *   a vector<X> or string branch b <p>.<T>.<b>..length[<p>.<T>..length]   Int_t, elements per row
 *                                  <p>.<T>.<b>..idx[<p>.<T>..length]      Int_t, start in the flat array
 *                                  <p>.<T>.<b>..totarraysize              Int_t, flat size
 *                                  <p>.<T>.<b>[<p>.<T>.<b>..totarraysize] the elements, X (Char_t for string)
 *   a vector<vector<X>> branch b   as vector<X> for the outer level (b..length/..idx/..totarraysize per row),
 *                                  then <p>.<T>.<b>.v..length / .v..idx [b..totarraysize],
 *                                  <p>.<T>.<b>.v..totarraysize, <p>.<T>.<b>.v[<p>.<T>.<b>.v..totarraysize]
 *
 * <p> is `prefix` (default "wc").  Values keep their exact ROOT type, so the
 * flat tree is a lossless re-layout of the row trees: row i of T, element j of
 * a vector branch b, is <p>.<T>.<b>[ <p>.<T>.<b>..idx[i] + j ].  An unsupported
 * branch type is an error (logged, and the event's flat tree is not written)
 * rather than a silent omission.
 *
 * The row trees are left in place unless drop_row_trees is true.
 */

#ifndef WIRECELLROOT_ROOTFLATTREEVISITOR
#define WIRECELLROOT_ROOTFLATTREEVISITOR

#include "WireCellClus/IEnsembleVisitor.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/Logging.h"

#include <string>
#include <vector>

namespace WireCell {
    namespace Root {

        class RootFlatTreeVisitor : public IConfigurable, public Clus::IEnsembleVisitor {
           public:
            RootFlatTreeVisitor();
            virtual ~RootFlatTreeVisitor();

            virtual void configure(const WireCell::Configuration& config);
            virtual Configuration default_configuration() const;
            virtual void visit(Clus::Facade::Ensemble& ensemble) const;

            /// Flatten every TTree of `filename` (except `tree_name` itself) into
            /// a one-entry tree `tree_name`.  Returns false (and writes nothing)
            /// on any unsupported branch type or I/O failure.  Public so tools
            /// and tests can apply it to an existing file.
            static bool flatten_file(const std::string& filename, const std::string& tree_name,
                                     const std::string& prefix, bool drop_row_trees,
                                     std::string& error);

           private:
            Log::logptr_t log;
            std::string m_output_filename{"tracking-pr.root"};
            std::string m_tree_name{"recTreeWireCell"};
            std::string m_prefix{"wc"};
            bool m_drop_row_trees{false};
        };

    }  // namespace Root
}  // namespace WireCell

#endif
