/** Tools to work with blobs
 */

#ifndef WIRECELLAUX_BLOBTOOLS
#define WIRECELLAUX_BLOBTOOLS

#include "WireCellIface/IBlobSet.h"

namespace WireCell::Aux {

    // Categorize a blob into okay or some type of malformation.
    struct BlobCategory {
        enum Status {okay, null, nostrips, nowidth, nocorners, novolume, numcats};
        Status stat;

        explicit BlobCategory(const IBlob::pointer& iblob);

        bool ok() const { return stat == okay; }
        static size_t size() {
            return (size_t) numcats;
        }
        size_t num() const {
            return (size_t) stat;
        }

        static size_t to_num(Status stat) {
            return static_cast<size_t>(stat);
        }
        static Status to_status(size_t index) {
            return static_cast<Status>(index);
        }
        static std::string to_str(size_t index)
        {
            return to_str(to_status(index));
        }

        static std::string to_str(Status stat)
        {
            switch (stat) {
                case okay: return "okay";
                case null: return "null";
                case nostrips: return "nostrips";
                case nowidth: return "nowidth";
                case nocorners: return "nocorners";
                case novolume: return "novolume";
                default: break;
            }
            return "undefined";
        }


        std::string str() const {
            return to_str(stat);
        }
    };

    // Dump a summary of a cluster graph.  If fingerprint, include more verbosity of channel-level info.
    std::string dumps(const IBlobSet::pointer& bs);

}
#endif
