#include "WireCellClus/MyFCN.h"
#include "WireCellUtil/Units.h"

#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>

#include <iostream>
#include <cmath>

using namespace WireCell;
using namespace WireCell::Clus::PR;

MyFCN::MyFCN(VertexPtr vtx, bool flag_vtx_constraint, double vtx_constraint_range, 
             double vertex_protect_dis, double vertex_protect_dis_short_track, double fit_dis)
    : vtx(vtx)
    , enforce_two_track_fit(false)
    , flag_vtx_constraint(flag_vtx_constraint)
    , vtx_constraint_range(vtx_constraint_range)
    , vertex_protect_dis(vertex_protect_dis)
    , vertex_protect_dis_short_track(vertex_protect_dis_short_track)
    , fit_dis(fit_dis)
{
    segments.clear();
    vec_points.clear();
}

MyFCN::~MyFCN()
{
}

void MyFCN::print_points()
{
    for (size_t i = 0; i != vec_points.size(); i++) {
        for (size_t j = 0; j != vec_points.at(i).size(); j++) {
            std::cout << i << " " << j << " " 
                      << vec_points.at(i).at(j).x() / units::cm << " " 
                      << vec_points.at(i).at(j).y() / units::cm << " " 
                      << vec_points.at(i).at(j).z() / units::cm << std::endl;
        }
    }
}

void MyFCN::AddSegment(SegmentPtr sg)
{
    // push in ...
    segments.push_back(sg);
    {
        std::vector<Facade::geo_point_t> pts;
        vec_points.push_back(pts);
    }

    Facade::geo_point_t center(0, 0, 0);
    double min_dis = 1e9;

    // Get fit points from segment
    const auto& fits = sg->fits();
    if (fits.empty()) {
        Facade::geo_point_t a(0, 0, 0);
        vec_PCA_dirs.push_back(std::make_tuple(a, a, a));
        vec_PCA_vals.push_back(std::make_tuple(0, 0, 0));
        vec_centers.push_back(a);
        return;
    }

    std::vector<Facade::geo_point_t> pts;
    for (const auto& fit : fits) {
        pts.push_back(fit.point);
    }
    double length = 0;
    if (pts.size() > 1) {
        auto front = pts.front();
        auto back = pts.back();
        length = std::sqrt(std::pow(front.x() - back.x(), 2) + 
                          std::pow(front.y() - back.y(), 2) + 
                          std::pow(front.z() - back.z(), 2));
    }

    Facade::geo_point_t vtx_pt = vtx->fit().point;
    
    for (size_t i = 0; i != pts.size(); i++) {
        double dis_to_vertex = std::sqrt(std::pow(pts.at(i).x() - vtx_pt.x(), 2) + 
                                        std::pow(pts.at(i).y() - vtx_pt.y(), 2) + 
                                        std::pow(pts.at(i).z() - vtx_pt.z(), 2));
        if (length > 3.0 * units::cm) {
            if (dis_to_vertex < vertex_protect_dis || dis_to_vertex > fit_dis) continue;
        } else {
            if (dis_to_vertex < vertex_protect_dis_short_track || dis_to_vertex > fit_dis) continue;
        }
        
        vec_points.back().push_back(pts.at(i));
        if (dis_to_vertex < min_dis) {
            center = pts.at(i);
            min_dis = dis_to_vertex;
        }
    }

    // calculate the PCA ...
    if (vec_points.back().size() > 1) {
        int nsum = vec_points.back().size();
        
        // Eigen vectors ...
        std::vector<Facade::geo_point_t> PCA_axis(3, Facade::geo_point_t(0, 0, 0));
        double PCA_values[3] = {0, 0, 0};
        
        Eigen::Matrix3d cov_matrix = Eigen::Matrix3d::Zero();
        
        for (size_t k = 0; k != vec_points.back().size(); k++) {
            double dx = vec_points.back().at(k).x() - center.x();
            double dy = vec_points.back().at(k).y() - center.y();
            double dz = vec_points.back().at(k).z() - center.z();
            
            cov_matrix(0, 0) += dx * dx;
            cov_matrix(0, 1) += dx * dy;
            cov_matrix(0, 2) += dx * dz;
            cov_matrix(1, 1) += dy * dy;
            cov_matrix(1, 2) += dy * dz;
            cov_matrix(2, 2) += dz * dz;
        }
        
        cov_matrix(1, 0) = cov_matrix(0, 1);
        cov_matrix(2, 0) = cov_matrix(0, 2);
        cov_matrix(2, 1) = cov_matrix(1, 2);
        
        cov_matrix /= nsum;
        
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen_solver(cov_matrix);
        Eigen::Vector3d eigen_values = eigen_solver.eigenvalues();
        Eigen::Matrix3d eigen_vectors = eigen_solver.eigenvectors();
        
        PCA_values[0] = eigen_values(0) + std::pow(0.15 * units::cm, 2);
        PCA_values[1] = eigen_values(1) + std::pow(0.15 * units::cm, 2);
        PCA_values[2] = eigen_values(2) + std::pow(0.15 * units::cm, 2);
        
        for (int i = 0; i != 3; i++) {
            double norm = std::sqrt(eigen_vectors(0, i) * eigen_vectors(0, i) + 
                                   eigen_vectors(1, i) * eigen_vectors(1, i) + 
                                   eigen_vectors(2, i) * eigen_vectors(2, i));
            PCA_axis[i] = Facade::geo_point_t(eigen_vectors(0, i) / norm, 
                                     eigen_vectors(1, i) / norm, 
                                     eigen_vectors(2, i) / norm);
        }
        
        vec_PCA_dirs.push_back(std::make_tuple(PCA_axis[0], PCA_axis[1], PCA_axis[2]));
        vec_PCA_vals.push_back(std::make_tuple(PCA_values[0], PCA_values[1], PCA_values[2]));
        vec_centers.push_back(center);
        
    } else {
        Facade::geo_point_t a(0, 0, 0);
        vec_PCA_dirs.push_back(std::make_tuple(a, a, a));
        vec_PCA_vals.push_back(std::make_tuple(0, 0, 0));
        if (vec_points.back().size() == 1) {
            vec_centers.push_back(vec_points.back().back());
        } else {
            vec_centers.push_back(a);
        }
    }
}

void MyFCN::update_fit_range(double tmp_vertex_protect_dis, double tmp_vertex_protect_dis_short_track, double tmp_fit_dis)
{
    vertex_protect_dis = tmp_vertex_protect_dis;
    vertex_protect_dis_short_track = tmp_vertex_protect_dis_short_track;
    fit_dis = tmp_fit_dis;

    std::vector<SegmentPtr> tmp_segments = segments;
    segments.clear();
    vec_points.clear();
    for (auto it = tmp_segments.begin(); it != tmp_segments.end(); it++) {
        AddSegment(*it);
    }
}

int MyFCN::get_fittable_tracks()
{
    int ncount = 0;
    for (size_t i = 0; i != vec_points.size(); i++) {
        if (vec_points.at(i).size() > 1) ncount++;
    }
    return ncount;
}

std::pair<SegmentPtr, int> MyFCN::get_seg_info(int i)
{
    if (i < (int)segments.size()) {
        return std::make_pair(segments.at(i), vec_points.at(i).size());
    }
    return std::make_pair(nullptr, 0);
}

std::pair<bool, WireCell::Clus::Facade::geo_point_t> MyFCN::FitVertex()
{
    Facade::geo_point_t fit_pos = vtx->fit().point;
    bool fit_flag = false;

    int ntracks = get_fittable_tracks();
    int npoints = 0;

    int n_large_angles = 0;
    for (size_t i = 0; i != vec_PCA_vals.size(); i++) {
        Eigen::Vector3d dir1(std::get<0>(vec_PCA_dirs.at(i)).x(), 
                            std::get<0>(vec_PCA_dirs.at(i)).y(), 
                            std::get<0>(vec_PCA_dirs.at(i)).z());
        for (size_t j = i + 1; j < vec_PCA_vals.size(); j++) {
            Eigen::Vector3d dir2(std::get<0>(vec_PCA_dirs.at(j)).x(), 
                                std::get<0>(vec_PCA_dirs.at(j)).y(), 
                                std::get<0>(vec_PCA_dirs.at(j)).z());
            double angle = std::acos(dir1.dot(dir2)) * 180.0 / M_PI;
            if (angle > 15) n_large_angles++;
        }
    }

    if ((ntracks > 2 && n_large_angles > 1) || (ntracks >= 2 && enforce_two_track_fit && n_large_angles >= 1)) {

        // start the fit ...
        Eigen::VectorXd temp_pos_3D_init(3), temp_pos_3D(3); // to be fitted
        temp_pos_3D_init(0) = fit_pos.x();
        temp_pos_3D_init(1) = fit_pos.y();
        temp_pos_3D_init(2) = fit_pos.z();

        Eigen::Vector3d b(0, 0, 0);
        Eigen::MatrixXd A = Eigen::MatrixXd::Zero(3, 3);

        for (size_t i = 0; i != vec_PCA_vals.size(); i++) {
            if (std::get<0>(vec_PCA_vals.at(i)) > 0) {
                npoints += vec_points.at(i).size();

                // fill the matrix ... first row, second column
                Eigen::MatrixXd R(3, 3);
                R(0, 0) = 0; R(0, 1) = 0; R(0, 2) = 0;
                double val1 = std::sqrt(std::get<0>(vec_PCA_vals.at(i)) / std::get<1>(vec_PCA_vals.at(i)));
                R(1, 0) = val1 * std::get<1>(vec_PCA_dirs.at(i)).x();
                R(1, 1) = val1 * std::get<1>(vec_PCA_dirs.at(i)).y();
                R(1, 2) = val1 * std::get<1>(vec_PCA_dirs.at(i)).z();
                val1 = std::sqrt(std::get<0>(vec_PCA_vals.at(i)) / std::get<2>(vec_PCA_vals.at(i)));
                R(2, 0) = val1 * std::get<2>(vec_PCA_dirs.at(i)).x();
                R(2, 1) = val1 * std::get<2>(vec_PCA_dirs.at(i)).y();
                R(2, 2) = val1 * std::get<2>(vec_PCA_dirs.at(i)).z();

                Eigen::Vector3d data(vec_centers.at(i).x(), vec_centers.at(i).y(), vec_centers.at(i).z());
                data = R * data;
                Eigen::MatrixXd RT = R.transpose();

                b += RT * data;
                A += RT * R;
            }
        }

        // add constraint ...
        if (flag_vtx_constraint) {
            Eigen::MatrixXd R = Eigen::MatrixXd::Zero(3, 3);
            R(0, 0) = 1.0 / vtx_constraint_range * std::sqrt(npoints);
            R(1, 1) = 1.0 / vtx_constraint_range * std::sqrt(npoints);
            R(2, 2) = 1.0 / vtx_constraint_range * std::sqrt(npoints);

            Eigen::Vector3d data(fit_pos.x() / vtx_constraint_range * std::sqrt(npoints), 
                                fit_pos.y() / vtx_constraint_range * std::sqrt(npoints), 
                                fit_pos.z() / vtx_constraint_range * std::sqrt(npoints));
            Eigen::MatrixXd RT = R.transpose();
            b += RT * data;
            A += RT * R;
        }

        Eigen::BiCGSTAB<Eigen::MatrixXd> solver;
        solver.compute(A);
        temp_pos_3D = solver.solveWithGuess(b, temp_pos_3D_init);

        if (!std::isnan(solver.error())) {
            fit_pos = Facade::geo_point_t(temp_pos_3D(0), temp_pos_3D(1), temp_pos_3D(2));
            fit_flag = true;
        } else {
            std::cout << "Cluster: ";
            if (vtx->cluster()) {
                std::cout << vtx->cluster()->get_cluster_id();
            } else {
                std::cout << "unknown";
            }
            std::cout << " Fit Vertex Failed!" << std::endl;
        }
    }

    return std::make_pair(fit_flag, fit_pos);
}

