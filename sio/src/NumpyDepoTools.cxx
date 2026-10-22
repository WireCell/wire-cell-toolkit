#include "NumpyDepoTools.h"
#include "WireCellUtil/String.h"
#include "WireCellUtil/Array.h"
#include "WireCellUtil/NumpyHelper.h"
#include "WireCellAux/SimpleDepo.h"


using namespace WireCell;

bool Sio::NumpyDepoTools::load(const std::string& fname, int count,
                               WireCell::IDepo::vector &depos)
{
    // match names used by NDS
    const std::string data_name = String::format("depo_data_%d", count);
    const std::string info_name = String::format("depo_info_%d", count);

    Array::array_xxf data;
    try {
        WireCell::Numpy::load2d(data, data_name, fname);
    }
    catch (std::runtime_error& err) {
        return false;
    }

    if (data.cols() != 7) {
        return false;
    }
    const size_t ndatas = data.rows();

    Array::array_xxi info;
    try {
        WireCell::Numpy::load2d(info, info_name, fname);
    }
    catch (std::runtime_error& err) {
        return false;
    }
    if (info.cols() != 4) {
        return false;
    }
    const size_t ninfos = info.rows();

    if (ndatas != ninfos) {
        return false;
    }
    const size_t ndepos = ndatas;

    // size_t npositive = 0;

    // Files written before 2026-10-21 carry child == own row index on
    // every gen>0 row (see Aux::fill); then the child is the previous row.
    size_t nprior_rows = 0, nself = 0;
    for (size_t ind=0; ind < ndepos; ++ind) {
        if (info(ind, 2) > 0) {
            ++nprior_rows;
            if ((size_t)info(ind, 3) == ind) ++nself;
        }
    }
    const bool legacy_child = nprior_rows > 0 && nself == nprior_rows;

    // All rows, active and prior, so that prior chains of any depth
    // can be linked (a gen 2 row's child is a gen 1 row).
    std::vector<std::shared_ptr<Aux::SimpleDepo>> all_sdepos;
    for (size_t ind=0; ind < ndepos; ++ind) {
        all_sdepos.push_back(std::make_shared<Aux::SimpleDepo>(
            data(ind, 0),        // t
            Point(data(ind, 2),  // x
                  data(ind, 3),  // y
                  data(ind, 4)), // z
            data(ind, 1),        // q
            nullptr,             // prior
            data(ind, 5),        // extent_long
            data(ind, 6),        // extent_tran
            info(ind, 0),        // id
            info(ind, 1)));      // pdg
    }

    for (size_t ind=0; ind < ndepos; ++ind) {
        const auto gen = info(ind, 2);
        if (gen == 0) {
            depos.push_back(all_sdepos[ind]);
            continue;
        }
        // this depo is a prior of row "other"
        size_t other = info(ind, 3);
        if (legacy_child) {
            other = ind - 1;
        }
        if (other >= all_sdepos.size() || other == ind) {
            return false;
        }
        all_sdepos[other]->set_prior(all_sdepos[ind]);
    }

    // (void)npositive; // silence unused warning

    return true;
}
