#include "WireCellAux/DepoTools.h"

using namespace WireCell;


// perhaps this should be exported?
struct DepoGenChild {
    IDepo::pointer depo;
    size_t gen, child;
};
using DepoGenChilds = std::vector<DepoGenChild>;


// Push depo and then, recursively, its prior chain.  The "child" of a
// prior row is the row index of the depo whose prior() it is, so that
// a reader can call child->set_prior(prior).  A gen==0 (active) row
// has child==0 by convention (unused).
//
// Defect fixed 2026-10-21 (wcp-porting-img pr_ml doc 19 sec 3.3): the
// recursion used to pass dp.size() taken AFTER the emplace, i.e. the
// index the prior itself was about to occupy, so every gen>0 row
// pointed at itself and DepoFileSource silently restored no priors.
static void push_depo(DepoGenChilds& dp,
                      IDepo::pointer depo,
                      size_t gen = 0, size_t childid = 0)
{
    dp.emplace_back(DepoGenChild{depo, gen, childid});
    const size_t myind = dp.size() - 1;
    auto prior = depo->prior();
    if (!prior) {
        return;
    }
    push_depo(dp, prior, gen + 1, myind);
}


static
DepoGenChilds flatten(const IDepo::vector& depos)
{
    DepoGenChilds ret;
    for (auto depo : depos) {
        push_depo(ret, depo);
    }
    return ret;
}


static
DepoGenChilds fake_flatten(const IDepo::vector& depos)
{
    DepoGenChilds ret;
    for (auto depo : depos) {
        ret.push_back(DepoGenChild{depo, 0,0});
    }
    return ret;
}


void Aux::fill(Array::array_xxf& data, Array::array_xxi& info,
               const IDepo::vector& your_depos, bool prior)
{
    DepoGenChilds fdepos;
    if (prior) {
        fdepos = flatten(your_depos);
    }
    else {
        fdepos = fake_flatten(your_depos);
    }

    const size_t nfdepos = fdepos.size();

    size_t nrows = fdepos.size();

    // time, charge, x, y, z, dlong, dtran
    const size_t ndata = 7;
    data = Array::array_xxf::Zero(nrows, ndata);

    // ID, pdg, gen, child
    const size_t ninfo = 4;
    info = Array::array_xxi::Zero(nrows, ninfo);

    for (size_t idepo = 0; idepo != nfdepos; ++idepo) {
        const auto& depogc = fdepos[idepo];
        const auto& depo = depogc.depo;
        data(idepo, 0) = depo->time();
        data(idepo, 1) = depo->charge();
        data(idepo, 2) = depo->pos().x();
        data(idepo, 3) = depo->pos().y();
        data(idepo, 4) = depo->pos().z();
        data(idepo, 5) = depo->extent_long();
        data(idepo, 6) = depo->extent_tran();
        info(idepo, 0) = depo->id();
        info(idepo, 1) = depo->pdg();
        info(idepo, 2) = depogc.gen;
        info(idepo, 3) = depogc.child;
    }

}

void Aux::fill_truth(Eigen::ArrayXXd& truth, const IDepo::vector& your_depos, bool prior)
{
    DepoGenChilds fdepos;
    if (prior) {
        fdepos = flatten(your_depos);
    }
    else {
        fdepos = fake_flatten(your_depos);
    }
    const size_t nrows = fdepos.size();

    // energy, orig_id, reserved
    const size_t ntruth = 3;
    truth = Eigen::ArrayXXd::Zero(nrows, ntruth);
    for (size_t idepo = 0; idepo != nrows; ++idepo) {
        truth(idepo, 0) = fdepos[idepo].depo->energy();
        // IDepo carries no pre-merge (G4 "orig") id; the column is
        // reserved for a Stage-A side table (doc 19 sec 3.3 item 3).
        truth(idepo, 1) = 0.0;
        truth(idepo, 2) = 0.0;
    }
}

std::vector<IDepo::vector>
Aux::depos_by_slice(const IDepo::vector& depos, const Binning& tbins,
                    double time_offset, double nsigma)
{
    std::vector<IDepo::vector> ret(tbins.nbins());
    for (const auto& depo : depos) {
        const double t = depo->time() + time_offset;
        const double sigma = depo->extent_long();
        const double tmin = t - sigma*nsigma;
        const double tmax = t + sigma*nsigma;
        int beg = std::max(tbins.bin(tmin),   0);
        const int end = std::min(tbins.bin(tmax)+1, tbins.nbins());
        while (beg < end) {
            ret[beg].push_back(depo);
            ++beg;
        }
    }
    return ret;
}

IDepo::vector Aux::sensitive(const IDepo::vector& depos, IAnodeFace::pointer face)
{
    IDepo::vector ret;
    auto bb = face->sensitive();
    if (bb.empty()) {
        return ret;
    }
    for (const auto& depo : depos) {
        if (bb.inside(depo->pos())) {
            ret.push_back(depo);
        }
    }
    return ret;
}
IDepo::vector Aux::sensitive(const IDepo::vector& depos, IAnodePlane::pointer anode)
{
    IDepo::vector ret;
    for (const auto& face : anode->faces()) {
        auto one = sensitive(depos, face);
        ret.insert(ret.end(), one.begin(), one.end());
    }
    return ret;
}
