#include <cmath>
#include <Eigen/Core>

//libigl includes
#include <igl/adjacency_list.h>
#include <igl/centroid.h>
#include <igl/bounding_box_diagonal.h>

#include "mcf_skeleton.h"

using namespace Eigen;
using ArrayXb = Eigen::Array<bool, Eigen::Dynamic,1>;

MCFSkeleton::MCFSkeleton(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const Eigen::MatrixXd& V_medial) 
    : V(V), F(F), MeanCurvatureFlow(V, F, V_medial) {
        
    setParameters(1.0, 20.0, 40.0, 0.002, 1e-6, 110);

    // Input vertices initially correspond to themselves
    corr = VectorXi::LinSpaced(V.rows(), 0, V.rows() - 1);
}

void MCFSkeleton::skeletonize(int max_iterations) { 
    Vector3d centroid;
    double volume;
    igl::centroid(V_temp, F_temp, centroid, volume);
    std::cout << "Initial Volume: " << volume << std::endl;

    iter = 0;
    while (volume > zero_TH && iter < max_iterations) {
        updateLaplacian();
        updateWeights();
        performMeshContraction();
        collapseShortEdges();
        splitBadTriangles();
        detectDegeneracies();

        igl::centroid(V_temp, F_temp, centroid, volume);
        std::cout << "Volume after" << " " << iter << " " << "iterations: " << volume << std::endl;
        iter++;
    }
}

// void MCFSkeleton::detectDegeneracies() {
//     std::vector<std::vector<int>> A;
//     igl::adjacency_list(F_temp, A);
//     for (int v = 0; v < V_temp.rows(); v++) {
//         int degenerate_count = 0;
//         for (int i = 0; i < A[v].size(); i++) {
//             int neighbor = A[v][i];
//             double mod_e = (V_temp.row(v) - V_temp.row(neighbor)).norm();
//             if (mod_e < 0.1 * edge_TH) {
//                 degenerate_count++;  
//             }
//             if (degenerate_count >= 2) {
//                 is_fixed(v) = true;
//                 if (is_split(v) && is_fixed(v)){
//                     std::cout << "Found clash" << std::endl;
//                 }
//                 continue;
//             }
//         }
//     }
// }

void MCFSkeleton::setParameters(double omega_L, double omega_H, double omega_M, double scale, double zero_TH, double alpha_TH) {
    this->omega_L = omega_L;
    this->omega_H = omega_H;
    this->omega_M = omega_M;
    this->edge_TH = scale * igl::bounding_box_diagonal(V);
    this->zero_TH = zero_TH;
    this->alpha_TH = alpha_TH / 180.0 * M_PI;
}