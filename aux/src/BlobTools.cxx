#include "WireCellAux/BlobTools.h"
#include "WireCellUtil/Units.h"

using namespace WireCell;

Aux::BlobCategory::BlobCategory(const IBlob::pointer& iblob)
{
    if (!iblob) {
        stat = Status::null;
        return;
    }

    const auto& blob = iblob->shape();
    const auto& strips = blob.strips();
    if (strips.empty()) {
        stat = Status::nostrips;
        return;
    }

    for (const auto& strip : strips) {
        if (strip.bounds.second == strip.bounds.first) {
            stat = Status::nowidth;
            return;
        }
    }
    const auto& corners = blob.corners();
    if (corners.empty()) {
        stat = Status::nocorners;
        return;
    }
    if (corners.size() < 3) {
        stat = Status::novolume;
        return;
    }
    stat = Status::okay;
}


#include <unordered_set>
std::string Aux::dumps(const IBlobSet::pointer& bs)
{
    std::stringstream ss;

    const size_t nbcats = BlobCategory::size();
    std::vector<size_t> blob_cats(nbcats, 0);
    std::unordered_set<ISlice::pointer> slices;

    for (auto iblob: bs->blobs()) {
        Aux::BlobCategory bcat(iblob);
        ++blob_cats[bcat.num()];
        slices.insert(iblob->slice());
    }

    ss << "Blobset ident:" << bs->ident()
       << " spans "  << slices.size() << " slices, main slice:";
    auto islice = bs->slice();
    if (islice) {
        ss << islice->ident() << " with " << islice->activity().size() << " activities, span="
           << islice->start()/units::us << "+" << islice->span()/units::us << " us";
        
    }
    else {
        ss << "(none)";
    }
    ss << " bcats:[";
    for (size_t icat=0; icat < nbcats; ++icat) {
        ss << " " << BlobCategory::to_str(icat) << "=" << blob_cats[icat];
    }
    ss << " ]";
    return ss.str();
}
