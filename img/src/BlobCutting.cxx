#include "WireCellImg/BlobCutting.h"
#include "WireCellAux/SimpleBlob.h"
#include "WireCellUtil/RayTiling.h"
#include "WireCellUtil/RayClustering.h"  // For overlap functions

#include "WireCellUtil/NamedFactory.h"

#include <iostream>
#include <algorithm>
#include <map>

WIRECELL_FACTORY(BlobCutting, WireCell::Img::BlobCutting,
                 WireCell::INamed,
                 WireCell::IFunctionNode<WireCell::IBlobSet, WireCell::IBlobSet>,
                 WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::RayGrid;

Img::BlobCutting::BlobCutting()
    : Aux::Logger("BlobCutting", "img")
{
}

Img::BlobCutting::~BlobCutting()
{
}

WireCell::Configuration Img::BlobCutting::default_configuration() const
{
    Configuration cfg;
    cfg["length_threshold"] = 20;    // Default threshold for strip length
    cfg["nudge"] = 0.01;            // Numerical precision parameter
    cfg["max_depth"] = 10;          // Maximum recursion depth
    return cfg;
}

void Img::BlobCutting::configure(const WireCell::Configuration& cfg)
{
    m_length_threshold = get(cfg, "length_threshold", 20);
    m_nudge = get(cfg, "nudge", 0.01);
    m_max_depth = get(cfg, "max_depth", 10);
}

// Helper function to check if a blob needs cutting
bool blob_needs_cutting(const RayGrid::Blob& blob, int length_threshold) {
    const auto& strips = blob.strips();
    for (const auto& strip : strips) {
        if (strip.layer >= 2 && strip.layer <= 4) { // U, V, W planes only
            int length = strip.bounds.second - strip.bounds.first;
            if (length > length_threshold) {
                return true;
            }
        }
    }
    return false;
}

// Helper function to check if a small blob is covered by (contained within) a big blob
// Uses the existing WireCell::RayGrid::surrounding() function for validation
bool is_blob_covered_by(const RayGrid::Blob& small_blob, const RayGrid::Blob& big_blob) {
    // const auto& small_strips = small_blob.strips();
    // const auto& big_strips = big_blob.strips();
    
    // // Create a map of layer -> strip for the big blob for easy lookup
    // std::map<int, RayGrid::Strip> big_strip_map;
    // for (const auto& strip : big_strips) {
    //     big_strip_map[strip.layer] = strip;
    // }
    
    // // Check if every strip in small blob is contained within corresponding strip in big blob
    // for (const auto& small_strip : small_strips) {
    //     auto big_it = big_strip_map.find(small_strip.layer);
    //     if (big_it == big_strip_map.end()) {
    //         // Small blob has a layer that big blob doesn't have - not covered
    //         return false;
    //     }
        
    //     const auto& big_strip = big_it->second;
        
    //     // Check if small strip bounds are within big strip bounds
    //     if (small_strip.bounds.first < big_strip.bounds.first || 
    //         small_strip.bounds.second > big_strip.bounds.second) {
    //         // Small strip extends beyond big strip - not covered
    //         return false;
    //     }
    // }
    
    // return true;
    // Additional validation using WireCell's surrounding() function
    // Convert to pointers for the existing function
    // Create blob vectors to get proper iterators for the surrounding() function
    RayGrid::blobs_t big_blobs = {big_blob};
    RayGrid::blobs_t small_blobs = {small_blob};
    
    // Get iterators (blobref_t) for the surrounding() function
    auto big_ref = big_blobs.begin();
    auto small_ref = small_blobs.begin();
    
    // Use WireCell's existing surrounding() function
    // surrounding(a, b) returns true if one blob completely surrounds the other
    // We want to check if big_blob surrounds small_blob
    bool is_surrounded = RayGrid::surrounding(big_ref, small_ref);
   
    return is_surrounded;



}

void print_blob(const RayGrid::Blob& blob) {
    std::cout << "Blob with " << blob.strips().size() << " strips and " 
              << blob.corners().size() << " corners:" << std::endl;
    for (const auto& strip : blob.strips()) {
        std::cout << "  Layer " << strip.layer 
                  << " bounds: [" << strip.bounds.first << ", " << strip.bounds.second << "]" << std::endl;
    }
    for (const auto& corner : blob.corners()) {
        std::cout << "  Corner: (L" << corner.first.layer << "@" << corner.first.grid 
                  << ", L" << corner.second.layer << "@" << corner.second.grid << ")" << std::endl;
    }   
}

// Function to split a blob into two parts along its longest strip
// The two parts together should cover exactly the same region as the original
std::vector<RayGrid::Blob> split_blob_once(
    const RayGrid::Coordinates& coords,
    const RayGrid::Blob& original_blob,
    double nudge)
{
    std::vector<RayGrid::Blob> result;
    
    const auto& strips = original_blob.strips();
    if (strips.size() < 3) {
        result.push_back(original_blob); // Can't split, return original
        return result;
    }
    
    // Find the longest strip among wire planes (layers 2, 3, 4)
    int longest_strip_index = -1;
    int max_length = 0;
    
    for (size_t i = 0; i < strips.size(); ++i) {
        const auto& strip = strips[i];
        if (strip.layer >= 2 && strip.layer <= 4) { // U, V, W planes only
            int length = strip.bounds.second - strip.bounds.first;
            if (length > max_length) {
                max_length = length;
                longest_strip_index = i;
            }
        }
    }
    
    if (longest_strip_index == -1 || max_length <= 2) {
        result.push_back(original_blob); // Can't split, return original
        return result;
    }
    
    const auto& longest_strip = strips[longest_strip_index];
    
    // Calculate the midpoint for division
    int mid_point = (longest_strip.bounds.first + longest_strip.bounds.second) / 2;
    
    // Create two strips that split the longest strip
    // IMPORTANT: These should be adjacent, not overlapping
    RayGrid::Strip first_half_strip{
        longest_strip.layer,
        {longest_strip.bounds.first, mid_point}
    };
    
    RayGrid::Strip second_half_strip{
        longest_strip.layer,
        {mid_point, longest_strip.bounds.second}
    };
    
    // Create first blob - use all original strips except replace the longest one
    RayGrid::Blob first_blob;
    for (int i = 0; i < (int)strips.size(); ++i) {
        if (i != longest_strip_index) {
            first_blob.add(coords, strips[i], nudge);
        }else{
            first_blob.add(coords, first_half_strip, nudge);
        }

    }
    
    // Create second blob - use all original strips except replace the longest one
    RayGrid::Blob second_blob;
    for (int i = 0; i < (int)strips.size(); ++i) {
        if (i != longest_strip_index) {
            second_blob.add(coords, strips[i], nudge);
        }else{
            second_blob.add(coords, second_half_strip, nudge);
        }
    }
    
    // Collect the new blobs for processing
    RayGrid::blobs_t new_blobs;
    if (!first_blob.corners().empty()) {
        new_blobs.push_back(first_blob);
    }
    if (!second_blob.corners().empty()) {
        new_blobs.push_back(second_blob);
    }
        std::cout << "ori blob:" << std::endl;
        print_blob(new_blobs[0]);
        print_blob(new_blobs[1]);
    
    // Apply WireCell's blob processing functions to ensure proper geometry
    // 1. Remove invalid blobs (those with empty corners or degenerate geometry)
    size_t dropped = RayGrid::drop_invalid(new_blobs);
    if (dropped > 0) {
        std::cout << "Dropped " << dropped << " invalid blobs after splitting" << std::endl;
    }
    

    // 2. Prune blob boundaries to actual geometric extents
    RayGrid::prune(coords, new_blobs, nudge);
    
    // 3. Drop invalid blobs again after pruning (pruning might create invalid geometry)
    dropped = RayGrid::drop_invalid(new_blobs);
    if (dropped > 0) {
        std::cout << "Dropped " << dropped << " invalid blobs after pruning" << std::endl;
    }
    
    // Now validate and add results - CHECK IF NEW BLOBS ARE COVERED BY ORIGINAL
    for (const auto& processed_blob : new_blobs) {
        std::cout << "Processed split blob:" << std::endl;
        print_blob(processed_blob);
        if (processed_blob.valid() && is_blob_covered_by(processed_blob, original_blob)) {
            result.push_back(processed_blob);
        } else {
            std::cout << "Processed split blob extends beyond original or is invalid - discarding" << std::endl;
        }
    }
    
    // If splitting failed, return original
    if (result.empty()) {
        result.push_back(original_blob);
    }
    
    return result;
}

// Recursive function to split blobs until all strips are below threshold
std::vector<RayGrid::Blob> split_blob_recursively(
    const RayGrid::Coordinates& coords,
    const RayGrid::Blob& blob,
    int length_threshold,
    double nudge,
    int max_depth)
{
    std::vector<RayGrid::Blob> result;
    
    // Base case: if blob doesn't need cutting or max depth reached
    if (!blob_needs_cutting(blob, length_threshold) || max_depth <= 0) {
        result.push_back(blob);
        return result;
    }
    
    // Split the blob once
    auto split_result = split_blob_once(coords, blob, nudge);
    
    if (split_result.size() <= 1) {
        // Splitting failed or produced only one blob, return it
        result.push_back(blob);
        return result;
    }
    
    // Recursively split each resulting blob
    for (const auto& sub_blob : split_result) {
        auto sub_results = split_blob_recursively(coords, sub_blob, length_threshold, nudge, max_depth - 1);
        result.insert(result.end(), sub_results.begin(), sub_results.end());
    }
    
    return result;
}

bool Img::BlobCutting::operator()(const input_pointer& in, output_pointer& out)
{
    // Handle end-of-stream
    if (!in) {
        out = nullptr;
        log->debug("EOS");
        return true;
    }

    // Create a new SimpleBlobSet for output
    auto sbs = std::make_shared<Aux::SimpleBlobSet>(in->ident(), in->slice());
    
    // Loop over all blobs in the input blob set
    const auto& input_blobs = in->blobs();
    // log->debug("BlobCutting: processing {} input blobs", input_blobs.size());
    
    int new_blob_id = 0; // Counter for new blob identifiers
    
    for (const auto& iblob : input_blobs) {
        const auto& shape = iblob->shape();
        
        if (blob_needs_cutting(shape, m_length_threshold)) {
            // Recursively split the blob until all strips are below threshold
            auto split_shapes = split_blob_recursively(
                iblob->face()->raygrid(), shape, m_length_threshold, m_nudge, m_max_depth);
            
            log->debug("Blob {} needed cutting, split into {} parts", 
                      iblob->ident(), split_shapes.size());
            
            // Create new IBlob objects for each split shape
            for (const auto& new_shape : split_shapes) {
                auto new_iblob = std::make_shared<Aux::SimpleBlob>(
                    new_blob_id++,  // Assign new unique identifier
                    iblob->value() / split_shapes.size(),  // Split the charge equally
                    iblob->uncertainty(),
                    new_shape,
                    iblob->slice(),
                    iblob->face()
                );
                sbs->m_blobs.push_back(new_iblob);
            }
        } else {
            // Blob doesn't need cutting, keep as-is
            sbs->m_blobs.push_back(iblob);
        }
    }
    
    out = sbs;
    // log->debug("BlobCutting: processed {} input blobs into {} output blobs", 
    //            input_blobs.size(), sbs->blobs().size());

    return true;
}