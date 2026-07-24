// SBND counterpart of wire-cell-uboone-magnify-tracking-convert.
//
// Fork by duplication (CLAUDE.md M10): the uBooNE app has other consumers
// (the qlport gate chain) and stays byte-for-byte untouched.  Differences
// from the uBooNE original, all of them SBND-driven:
//
//   1. NO hard-coded uBooNE x transform.  The uBooNE app rewrites truth x as
//      (x+0.6)/1.098*1.1009999-0.1101 (its lines 148/183) to undo a
//      simulation-vs-imaging drift-speed mismatch plus a Y/U plane offset.
//      Those constants are uBooNE's.  Here the truth x is passed through
//      unchanged by default; -s<scale> -c<offset> apply x -> x*scale + offset
//      if an SBND-specific correction is ever needed.  This is what makes MC
//      mode (-f1) usable for SBND truth pairing at all.
//   2. Carries the STM dump's extra per-point branches (pass, status, rr)
//      into T_rec as nested vectors, so a Magnify-tracking file made from
//      SbndMagnifyTrackingVisitor output keeps the per-pass provenance.
//   3. Clones the STM decision trees T_stm_pass / T_stm_eval into the output
//      when present, keeping the Magnify file self-contained.
//
// Input is the tracking ROOT file written by SbndMagnifyTrackingVisitor
// (T_rec_charge / T_proj_data / T_bad_ch / Trun, plus T_stm_pass /
// T_stm_eval); output is the five-tree format the Magnify-tracking GUI reads
// (T_rec / T_proj_data / T_proj / T_bad_ch / optional T_true).
//
// Block ids: SbndMagnifyTrackingVisitor writes ndf = cluster_id*10 + pass, so
// forward and backward fits of one cluster arrive as separate tracks.  As in
// the uBooNE app, a new track starts whenever ndf changes, which assumes the
// writer emits each block contiguously (it does).

#include "WireCellUtil/NFKDVec.h"
#include "WireCellUtil/Point.h"
#include "WireCellUtil/Units.h"

#include <iostream>
#include <vector>
#include <map>

#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include "TString.h"
#include "TVector3.h"

using namespace WireCell;

// Helper: add a PointVector to an NFKDVec tree
static void add_points(NFKDVec::Tree<double>& tree, const PointVector& ps) {
    if (ps.empty()) return;
    std::vector<std::vector<double>> coords(3);
    for (int i = 0; i < 3; ++i) coords[i].reserve(ps.size());
    for (const auto& p : ps) {
        coords[0].push_back(p.x());
        coords[1].push_back(p.y());
        coords[2].push_back(p.z());
    }
    tree.append(coords);
}

// Helper: find the closest point in the tree, returns (distance, closest_point)
static std::pair<double, Point> get_closest_point(const NFKDVec::Tree<double>& tree, const Point& p) {
    auto results = tree.knn(1, p);
    if (results.empty()) return {0, Point(0, 0, 0)};
    size_t idx = results[0].first;
    double dist = std::sqrt(results[0].second);  // L2Simple returns squared distance
    return {dist, tree.point3d(idx)};
}

int main(int argc, char* argv[])
{
    if (argc < 2) {
        std::cerr << "usage: wire-cell-sbnd-magnify-tracking-convert "
                  << "-b[reco.root] -t[reco_treename] -o[out.root] -f[1:MC,2:data] "
                  << "[-a truth.root -n truth_treename] [-s x_scale] [-c x_offset]"
                  << std::endl;
        std::cerr << "  MC mode (-f1) pairs each fitted point with the nearest truth point;"
                  << " truth x is passed through unchanged unless -s/-c are given." << std::endl;
        return 1;
    }

    TString reco_filename = "tracking-stm.root";
    TString reco_treename = "T_rec_charge";
    TString proj_treename = "T_proj_data";

    TString truth_filename = "mcs-tracks.root";
    TString truth_treename = "T";
    TString out_filename = "track_com.root";

    int file_type = 2;             // 1 = MC (truth pairing), 2 = data
    double x_scale = 1.0;          // truth x -> x*x_scale + x_offset (identity by default)
    double x_offset = 0.0;

    for (Int_t i = 1; i != argc; i++) {
        switch (argv[i][1]) {
        case 'b': reco_filename  = &argv[i][2]; break;
        case 't': reco_treename  = &argv[i][2]; break;
        case 'a': truth_filename = &argv[i][2]; break;
        case 'n': truth_treename = &argv[i][2]; break;
        case 'o': out_filename   = &argv[i][2]; break;
        case 'f': file_type      = atoi(&argv[i][2]); break;
        case 's': x_scale        = atof(&argv[i][2]); break;
        case 'c': x_offset       = atof(&argv[i][2]); break;
        }
    }

    if (file_type == 1) {
        std::cout << truth_filename << " " << reco_filename << " "
                  << truth_treename << " " << reco_treename
                  << "  (truth x -> x*" << x_scale << " + " << x_offset << ")" << std::endl;
    }
    else {
        std::cout << reco_filename << " " << reco_treename << std::endl;
    }

    TFile* file1 = new TFile(reco_filename);
    if (file1->IsZombie()) {
        std::cerr << "ERROR: cannot open " << reco_filename << std::endl;
        return 1;
    }
    TTree* T_bad_ch = (TTree*) file1->Get("T_bad_ch");
    TTree* T_stm_pass = (TTree*) file1->Get("T_stm_pass");
    TTree* T_stm_eval = (TTree*) file1->Get("T_stm_eval");

    Double_t dQdx_scale = 1., dQdx_offset = 0.;
    TTree* Trun = (TTree*) file1->Get("Trun");
    if (Trun != 0 && Trun->GetBranch("dQdx_scale")) {
        Trun->SetBranchAddress("dQdx_scale", &dQdx_scale);
        Trun->SetBranchAddress("dQdx_offset", &dQdx_offset);
        Trun->GetEntry(0);
    }

    TFile* file = new TFile(out_filename, "RECREATE");
    // Clones attach to the output file and are written once by file->Write()
    // below.  (The uBooNE app Write()s them here as well, which leaves two
    // cycles of each tree in its output; harmless but avoided here.)
    if (T_bad_ch != 0) T_bad_ch->CloneTree();
    // STM decision trees ride along so the Magnify file explains its own tracks.
    if (T_stm_pass != 0) T_stm_pass->CloneTree();
    if (T_stm_eval != 0) T_stm_eval->CloneTree();

    std::vector<std::vector<double> >* vx = new std::vector<std::vector<double> >;
    std::vector<std::vector<double> >* vy = new std::vector<std::vector<double> >;
    std::vector<std::vector<double> >* vz = new std::vector<std::vector<double> >;
    std::vector<std::vector<double> >* vQ = new std::vector<std::vector<double> >;
    std::vector<int>* vN = new std::vector<int>;

    std::vector<double>* x = new std::vector<double>;
    std::vector<double>* y = new std::vector<double>;
    std::vector<double>* z = new std::vector<double>;
    std::vector<double>* Q = new std::vector<double>;
    Int_t N;

    if (file_type == 1) {  // MC: read the truth track and republish it as T_true
        TChain* T_true = new TChain(truth_treename, truth_treename);
        Int_t nfiles = T_true->Add(truth_filename);

        if (nfiles > 0 && T_true->GetEntries() > 0) {
            TTree* t2 = new TTree("T_true", "T_true");
            T_true->SetBranchAddress("N", &N);
            T_true->SetBranchAddress("x", &x);
            T_true->SetBranchAddress("y", &y);
            T_true->SetBranchAddress("z", &z);
            T_true->SetBranchAddress("Q", &Q);
            T_true->GetEntry(0);
            vN->push_back(N);
            vx->push_back(*x);
            vy->push_back(*y);
            vz->push_back(*z);
            vQ->push_back(*Q);

            // Identity unless -s/-c were given.  No uBooNE constants here.
            for (size_t i = 0; i != vx->at(0).size(); i++) {
                vx->at(0).at(i) = vx->at(0).at(i) * x_scale + x_offset;
            }

            t2->Branch("N", &vN);
            t2->Branch("x", &vx);
            t2->Branch("y", &vy);
            t2->Branch("z", &vz);
            t2->Branch("Q", &vQ);
            t2->Fill();
        }
        else {
            std::cerr << "Warning: truth file " << truth_filename
                      << " does not exist or is empty.  Skipping truth tree." << std::endl;
        }
        delete T_true;
    }

    NFKDVec::Tree<double> pcloud(3), pcloud1(3);
    PointVector ps;

    if (file_type == 1 && !x->empty()) {
        for (size_t i = 0; i != x->size(); i++) {
            x->at(i) = x->at(i) * x_scale + x_offset;
            ps.push_back(Point(x->at(i) * units::cm, y->at(i) * units::cm, z->at(i) * units::cm));
        }
        add_points(pcloud, ps);
    }
    ps.clear();

    TChain* T_rec = new TChain(reco_treename, reco_treename);
    TChain* T_proj_data = new TChain(proj_treename, proj_treename);
    TChain* T_proj = new TChain("T_proj", "T_proj");
    T_rec->Add(reco_filename);
    T_proj_data->Add(reco_filename);
    T_proj->Add(reco_filename);

    Double_t x1, y1, z1;
    Double_t dQ1, dx1, ndf;
    Double_t pu, pv, pw, pt;
    T_rec->SetBranchAddress("x", &x1);
    T_rec->SetBranchAddress("y", &y1);
    T_rec->SetBranchAddress("z", &z1);
    T_rec->SetBranchAddress("q", &dQ1);
    T_rec->SetBranchAddress("nq", &dx1);
    T_rec->SetBranchAddress("ndf", &ndf);
    T_rec->SetBranchAddress("pu", &pu);
    T_rec->SetBranchAddress("pv", &pv);
    T_rec->SetBranchAddress("pw", &pw);
    T_rec->SetBranchAddress("pt", &pt);

    Double_t reduced_chi2 = 0;
    const bool has_chi2 = T_rec->GetBranch("reduced_chi2");
    if (has_chi2) T_rec->SetBranchAddress("reduced_chi2", &reduced_chi2);

    // STM extras written by SbndMagnifyTrackingVisitor.
    Double_t rr = 0;
    Int_t stm_pass = 0, stm_status = -1;
    const bool has_rr     = T_rec->GetBranch("rr");
    const bool has_pass   = T_rec->GetBranch("pass");
    const bool has_status = T_rec->GetBranch("status");
    if (has_rr)     T_rec->SetBranchAddress("rr", &rr);
    if (has_pass)   T_rec->SetBranchAddress("pass", &stm_pass);
    if (has_status) T_rec->SetBranchAddress("status", &stm_status);

    TTree* t1 = new TTree("T_rec", "T_rec");
    t1->SetDirectory(file);

    std::vector<double>* max_dis = new std::vector<double>;
    std::vector<double>* beg_dis = new std::vector<double>;
    std::vector<double>* end_dis = new std::vector<double>;
    std::vector<double>* total_dis2 = new std::vector<double>;
    std::vector<double>* total_L = new std::vector<double>;
    std::vector<int>* Npoints = new std::vector<int>;
    std::vector<double>* total_dtheta = new std::vector<double>;
    std::vector<double>* max_dtheta = new std::vector<double>;
    if (file_type == 1) {
        t1->Branch("stat_max_dis", &max_dis);
        t1->Branch("stat_beg_dis", &beg_dis);
        t1->Branch("stat_end_dis", &end_dis);
        t1->Branch("stat_total_dis2", &total_dis2);
        t1->Branch("stat_total_dtheta", &total_dtheta);
        t1->Branch("stat_max_dtheta", &max_dtheta);
    }
    t1->Branch("stat_total_L", &total_L);
    t1->Branch("stat_N", &Npoints);

    std::vector<std::vector<double> >* x2 = new std::vector<std::vector<double> >;
    std::vector<std::vector<double> >* y2 = new std::vector<std::vector<double> >;
    std::vector<std::vector<double> >* z2 = new std::vector<std::vector<double> >;
    std::vector<std::vector<double> >* dQ_rec = new std::vector<std::vector<double> >;
    std::vector<std::vector<double> >* dQ_tru = new std::vector<std::vector<double> >;
    std::vector<std::vector<double> >* dx = new std::vector<std::vector<double> >;
    std::vector<int>* cluster_id = new std::vector<int>;
    std::vector<std::vector<double> >* rec_pu = new std::vector<std::vector<double> >;
    std::vector<std::vector<double> >* rec_pv = new std::vector<std::vector<double> >;
    std::vector<std::vector<double> >* rec_pw = new std::vector<std::vector<double> >;
    std::vector<std::vector<double> >* rec_pt = new std::vector<std::vector<double> >;

    std::vector<std::vector<double> >* x2_pair = new std::vector<std::vector<double> >;
    std::vector<std::vector<double> >* y2_pair = new std::vector<std::vector<double> >;
    std::vector<std::vector<double> >* z2_pair = new std::vector<std::vector<double> >;

    std::vector<std::vector<double> >* dis = new std::vector<std::vector<double> >;
    std::vector<std::vector<double> >* L = new std::vector<std::vector<double> >;
    std::vector<std::vector<double> >* dtheta = new std::vector<std::vector<double> >;

    std::vector<std::vector<double> >* Vreduced_chi2 = new std::vector<std::vector<double> >;
    std::vector<std::vector<double> >* Vrr = new std::vector<std::vector<double> >;
    std::vector<std::vector<int> >* Vpass = new std::vector<std::vector<int> >;
    std::vector<std::vector<int> >* Vstatus = new std::vector<std::vector<int> >;

    t1->Branch("rec_x", &x2);
    t1->Branch("rec_y", &y2);
    t1->Branch("rec_z", &z2);
    t1->Branch("rec_dQ", &dQ_rec);
    t1->Branch("rec_dx", &dx);
    t1->Branch("rec_L", &L);
    t1->Branch("rec_cluster_id", &cluster_id);
    t1->Branch("rec_u", &rec_pu);
    t1->Branch("rec_v", &rec_pv);
    t1->Branch("rec_w", &rec_pw);
    t1->Branch("rec_t", &rec_pt);
    if (has_chi2)   t1->Branch("reduced_chi2", &Vreduced_chi2);
    if (has_rr)     t1->Branch("rec_rr", &Vrr);
    if (has_pass)   t1->Branch("stm_pass", &Vpass);
    if (has_status) t1->Branch("stm_status", &Vstatus);

    if (file_type == 1) {
        t1->Branch("true_dQ", &dQ_tru);
        t1->Branch("true_x", &x2_pair);
        t1->Branch("true_y", &y2_pair);
        t1->Branch("true_z", &z2_pair);
        t1->Branch("com_dis", &dis);
        t1->Branch("com_dtheta", &dtheta);
    }

    double prev_x1{0}, prev_y1{0}, prev_z1{0};
    int prev_cluster_id = -1;

    std::map<std::tuple<int, int, int>, std::pair<int, int> > map_point_index;

    for (int i = 0; i != T_rec->GetEntries(); i++) {
        T_rec->GetEntry(i);
        if (std::round(ndf) != prev_cluster_id) {   // new block = new track
            Npoints->push_back(0);
            total_L->push_back(0);

            dQ_tru->push_back(std::vector<double>());
            dis->push_back(std::vector<double>());
            x2_pair->push_back(std::vector<double>());
            y2_pair->push_back(std::vector<double>());
            z2_pair->push_back(std::vector<double>());

            x2->push_back(std::vector<double>());
            y2->push_back(std::vector<double>());
            z2->push_back(std::vector<double>());
            rec_pu->push_back(std::vector<double>());
            rec_pv->push_back(std::vector<double>());
            rec_pw->push_back(std::vector<double>());
            rec_pt->push_back(std::vector<double>());
            cluster_id->push_back(std::round(ndf));
            dQ_rec->push_back(std::vector<double>());
            dx->push_back(std::vector<double>());
            L->push_back(std::vector<double>());

            if (has_chi2)   Vreduced_chi2->push_back(std::vector<double>());
            if (has_rr)     Vrr->push_back(std::vector<double>());
            if (has_pass)   Vpass->push_back(std::vector<int>());
            if (has_status) Vstatus->push_back(std::vector<int>());

            max_dis->push_back(0);
            total_dis2->push_back(0);
        }
        Npoints->back()++;

        if (Npoints->back() != 1) {
            total_L->back() += sqrt(pow(x1 - prev_x1, 2) + pow(y1 - prev_y1, 2) + pow(z1 - prev_z1, 2));
        }

        Point p(x1 * units::cm, y1 * units::cm, z1 * units::cm);
        ps.push_back(p);

        if (file_type == 1 && !x->empty()) {
            std::pair<double, Point> point_pair = get_closest_point(pcloud, p);
            map_point_index[std::make_tuple(p.x() / (0.01 * units::mm),
                                            p.y() / (0.01 * units::mm),
                                            p.z() / (0.01 * units::mm))] =
                std::make_pair(x2->size() - 1, x2->back().size());
            dQ_tru->back().push_back(0);
            dis->back().push_back(point_pair.first / units::cm);

            x2_pair->back().push_back(point_pair.second.x() / units::cm);
            y2_pair->back().push_back(point_pair.second.y() / units::cm);
            z2_pair->back().push_back(point_pair.second.z() / units::cm);

            if (max_dis->back() < point_pair.first / units::cm)
                max_dis->back() = point_pair.first / units::cm;
            total_dis2->back() += pow(point_pair.first / units::cm, 2);

            const double dis1 = pow(x1 - x->front(), 2) + pow(y1 - y->front(), 2) + pow(z1 - z->front(), 2);
            const double dis2 = pow(x1 - x->back(), 2) + pow(y1 - y->back(), 2) + pow(z1 - z->back(), 2);
            if (Npoints->back() == 1) {
                beg_dis->push_back(sqrt(std::min(dis1, dis2)));
                end_dis->push_back(0);
            }
            end_dis->back() = sqrt(std::min(dis1, dis2));
        }

        x2->back().push_back(x1);
        y2->back().push_back(y1);
        z2->back().push_back(z1);
        rec_pu->back().push_back(pu);
        rec_pv->back().push_back(pv);
        rec_pw->back().push_back(pw);
        rec_pt->back().push_back(pt);
        dQ_rec->back().push_back((dQ1 - dQdx_offset) / dQdx_scale);
        dx->back().push_back(dx1);
        L->back().push_back(total_L->back());

        if (has_chi2)   Vreduced_chi2->back().push_back(reduced_chi2);
        if (has_rr)     Vrr->back().push_back(rr);
        if (has_pass)   Vpass->back().push_back(stm_pass);
        if (has_status) Vstatus->back().push_back(stm_status);

        prev_x1 = x1;
        prev_y1 = y1;
        prev_z1 = z1;
        prev_cluster_id = std::round(ndf);
    }

    add_points(pcloud1, ps);

    if (file_type == 1 && !x->empty()) {
        // Accumulate truth charge onto the fitted point each truth point is
        // closest to, then per-track angular comparison.
        for (size_t i = 0; i != x->size(); i++) {
            Point p(x->at(i) * units::cm, y->at(i) * units::cm, z->at(i) * units::cm);
            std::pair<double, Point> point_pair = get_closest_point(pcloud1, p);
            auto key = std::make_tuple(int(point_pair.second.x() / (0.01 * units::mm)),
                                       int(point_pair.second.y() / (0.01 * units::mm)),
                                       int(point_pair.second.z() / (0.01 * units::mm)));
            const int index = map_point_index[key].second;
            const int index1 = map_point_index[key].first;
            if (static_cast<size_t>(index1) < dQ_tru->size() &&
                static_cast<size_t>(index) < dQ_tru->at(index1).size())
                dQ_tru->at(index1).at(index) += Q->at(i);
        }

        for (size_t k = 0; k != x2->size(); k++) {
            if (x2->at(k).size() > 1) {
                dtheta->push_back(std::vector<double>());
                max_dtheta->push_back(0);
                total_dtheta->push_back(0);

                const size_t n = x2->at(k).size();
                for (size_t i = 0; i != n; i++) {
                    if (i == 0) {
                        TVector3 dir1(x2->at(k).at(1) - x2->at(k).at(0),
                                      y2->at(k).at(1) - y2->at(k).at(0),
                                      z2->at(k).at(1) - z2->at(k).at(0));
                        TVector3 dir2(x2_pair->at(k).at(1) - x2_pair->at(k).at(0),
                                      y2_pair->at(k).at(1) - y2_pair->at(k).at(0),
                                      z2_pair->at(k).at(1) - z2_pair->at(k).at(0));
                        dtheta->at(k).push_back(dir1.Angle(dir2));
                    }
                    else if (i == n - 1) {
                        TVector3 dir1(x2->at(k).at(n - 1) - x2->at(k).at(n - 2),
                                      y2->at(k).at(n - 1) - y2->at(k).at(n - 2),
                                      z2->at(k).at(n - 1) - z2->at(k).at(n - 2));
                        TVector3 dir2(x2_pair->at(k).at(n - 1) - x2_pair->at(k).at(n - 2),
                                      y2_pair->at(k).at(n - 1) - y2_pair->at(k).at(n - 2),
                                      z2_pair->at(k).at(n - 1) - z2_pair->at(k).at(n - 2));
                        dtheta->at(k).push_back(dir1.Angle(dir2));
                    }
                    else {
                        TVector3 dir1(x2->at(k).at(i + 1) - x2->at(k).at(i),
                                      y2->at(k).at(i + 1) - y2->at(k).at(i),
                                      z2->at(k).at(i + 1) - z2->at(k).at(i));
                        TVector3 dir2(x2_pair->at(k).at(i + 1) - x2_pair->at(k).at(i),
                                      y2_pair->at(k).at(i + 1) - y2_pair->at(k).at(i),
                                      z2_pair->at(k).at(i + 1) - z2_pair->at(k).at(i));
                        TVector3 dir3(x2->at(k).at(i - 1) - x2->at(k).at(i),
                                      y2->at(k).at(i - 1) - y2->at(k).at(i),
                                      z2->at(k).at(i - 1) - z2->at(k).at(i));
                        TVector3 dir4(x2_pair->at(k).at(i - 1) - x2_pair->at(k).at(i),
                                      y2_pair->at(k).at(i - 1) - y2_pair->at(k).at(i),
                                      z2_pair->at(k).at(i - 1) - z2_pair->at(k).at(i));
                        dtheta->at(k).push_back((dir1.Angle(dir2) + dir3.Angle(dir4)) / 2.);
                    }
                    if (dtheta->at(k).back() > max_dtheta->at(k))
                        max_dtheta->at(k) = dtheta->at(k).back();
                    total_dtheta->at(k) += dtheta->at(k).back();
                }
            }
            else {
                dtheta->push_back(std::vector<double>(1, 0.));
                max_dtheta->push_back(0);
                total_dtheta->push_back(0);
            }
        }
    }
    t1->Fill();

    T_proj_data->CloneTree(-1, "fast");
    T_proj->CloneTree(-1, "fast");

    file->Write();
    file->Close();

    std::cout << "wrote " << out_filename << ": " << cluster_id->size()
              << " track block(s), " << T_rec->GetEntries() << " fitted point(s)" << std::endl;
    return 0;
}
