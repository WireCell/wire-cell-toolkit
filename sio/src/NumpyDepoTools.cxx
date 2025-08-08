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

    std::vector<Aux::SimpleDepo*> sdepos;
    for (size_t ind=0; ind < ndepos; ++ind) {

        // {
        //     auto q = data(ind, 1);
        //     if (q > 0) {
        //         ++npositive;
        //     }
        // }

        auto sdepo = new Aux::SimpleDepo(
            data(ind, 0),        // t
            Point(data(ind, 2),  // x
                  data(ind, 3),  // y
                  data(ind, 4)), // z
            data(ind, 1),        // q
            nullptr,             // prior
            data(ind, 5),        // extent_long
            data(ind, 6),        // extent_tran
            info(ind, 0),        // id
            info(ind, 1));       // pdg

        const auto gen = info(ind, 2);
        if (gen > 0) {
            // this depo is a prior
            const size_t other = info(ind, 3);
            if (other >= sdepos.size()) {
                return false;
            }
        
            auto idepo = IDepo::pointer(sdepo);
            sdepo = nullptr;
            sdepos[other]->set_prior(idepo);
        }
        // We save both gen=0 and gen>0 as nullptrs to preserve indexing
        sdepos.push_back(sdepo);
    }
    
    for (auto sdepo: sdepos) {
        if (sdepo) {
            auto idepo = IDepo::pointer(sdepo);
            sdepo = nullptr;
            depos.push_back(idepo);
        }
    }

    return true;
}
