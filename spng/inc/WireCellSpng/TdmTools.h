#ifndef WIRECELL_SPNG_TDMTOOLS
#define WIRECELL_SPNG_TDMTOOLS

#include "WireCellSpng/ITorchTensor.h"
#include "WireCellUtil/Configuration.h"
#include <string>
#include <vector>

namespace WireCell::SPNG::TDM {

    /// Handle datapath and provenance metadata parameters.
    ///
    /// @param md A metadata object for an ITorchTensor under construction.
    /// @param from A metadata object for a given ITorchTensor from which the under-construction tensor is derived.
    /// @param datapath_format A string with "{fmt} style {params}" to form the "datapath" entry.
    /// @return A metadata object fit for including in the under-construction ITorchTensor.
    ///
    /// The output metadata object will include a "derived_from" attribute of
    /// value string holding the "datapath" attribute of from, if provided.
    Configuration derive_metadata(Configuration md,
                                  Configuration from,
                                  const std::string& datapath_format = "");


    /// An multi-derived metadata object.
    ///
    /// This works as the scalar version above except that it will merge
    /// metadata objects from multiple given tensors.  The merge conflict policy
    /// is last one wins.
    ///
    /// The "derived_from" will be an array of datapath strings.
    Configuration derive_metadata(Configuration md,
                                  const std::vector<Configuration>& from,
                                  const std::string& datapath_format = "");


    /**
       @brief Determine if an object matches.

       @param match A match object
       @param target An object to apply match algorithm.

       For target to match the following must be true.

       - The target object provides all the attribute key names provided by the
       match object.

       - Every target attribute value in this set is comparable to the match
       attribute value.

       Attribute values compare if any one of these are true:

       - The target and match attribute values are of the same type and equal.

       - The match attribute is a string and matches via regex to the target
         value after that is coerced to string type.

       Note, this matching is policy free.  There are no "special" attributes
       nor "special" matching.  
     */
    bool match_object(const Configuration& match, const Configuration& object);

    /**
       @brief Select tensor via object matching.

       @param tensors A vector of shared pointer to ITorchTensor from which to select.
       @param matches A JSON object or array of objects providing the match specification.
       @return A vector of shared pointer to ITorchTensor holding selection.

       The match is applied to each tensor metadata object via match_object()
       and all that match are returned.

       The selected tenors output order is determined by the order of matches
       (if it is an array) and the order of the tensors.  The order is
       tensor-major.  That is, each match object may produce a contiguous block
       of tensors in the output order.

       Note, as match_object() is policy free the user should take care to make
       a match unique to the desired purpose.  For example, if the goal is to
       select tensors of a particular TDM datatype, include a "datatype"
       attribute in the match object.

       Note, this selection does not scale well over large sets of tensors.
       Partitioning a large set given prior knowledge of the matches can help.
       See by_datatype() as one possibly useful partitioning.

     */
    ITorchTensor::vector select_tensors(const ITorchTensor::vector& tensors,
                                        const Configuration& matches);

    /**
       @brief Partition tensors by their datatype.

       @param tenors A set of tensors.
       @return A map from datatype to a set of tensors with that datatype.
       
       Any tensor lacking datatype will be placed in the map with an empty string key.
     */
    std::map<std::string, ITorchTensor::vector> by_datatype(const ITorchTensor::vector& tensors);

}

#endif
