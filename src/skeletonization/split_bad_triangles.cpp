#include <algorithm> 
#include <cassert>
#include <cmath>
#include <functional>
#include <Eigen/Core> 

// libigl includes                   
#include <igl/internal_angles.h>           
#include <igl/edge_lengths.h>             
#include <igl/edge_flaps.h>                
#include <igl/per_vertex_normals.h> 
#include <igl/unique_edge_map.h> 
#include <igl/unique.h>  
#include <igl/dot_row.h>    

#include "mcf_skeleton.h" // Class declaration and utility functions

using namespace Eigen;
using ArrayXb = Eigen::Array<bool, Eigen::Dynamic, 1>; // Alias for a dynamic-sized array of booleans

void MCFSkeleton::splitBadTriangles() {
    auto bad_triangles = [this] (const MatrixXd& V, const MatrixXi& F, const MatrixXi& /*E*/, 
            const MatrixXi& /*uE*/, const VectorXi& /*EMAP*/, const std::vector<std::vector<int>>& uE2E, 
            std::priority_queue<std::tuple<double,int,int>> &Q, int uei) -> bool {

            int num_faces = F.rows();
            int f0 = uE2E[uei][0] % num_faces;
            int f1 = uE2E[uei][1] % num_faces;
            int c0 = uE2E[uei][0] / num_faces;
            int c1 = uE2E[uei][1] / num_faces;
            int e0 = F(f0, (c0+1)%3); // first vertex of split edge
            int e1 = F(f0, (c0+2)%3); // second vertex of split edge
            int v0 = F(f0, c0);
            int v1 = F(f1, c1);

            if (is_fixed(e0) && is_fixed(e1)) {
                return false; // Do not split fixed edges
            } 

            double alpha0 = utils::compute_internal_angle(V, e0, e1, v0);
            double alpha1 = utils::compute_internal_angle(V, e0, e1, v1);
            if (alpha0 < alpha_TH && alpha1 < alpha_TH) {
                return false;
            }
            
            // Using index 2 entry to log projected vertex
            if ((V.row(e0) - V.row(v0)).norm() < edge_TH || 
                (V.row(e1) - V.row(v0)).norm() < edge_TH || 
                (V.row(e0) - V.row(v1)).norm() < edge_TH || 
                (V.row(e1) - V.row(v1)).norm() < edge_TH || 
                (V.row(e0) - V.row(e1)).norm() < edge_TH) {
                return false;  // Do not split degenerate triangles
            }

            std::tuple<double, int, int> bad_corner;
            bad_corner = (alpha0 > alpha1) ? std::make_tuple(alpha0, uei, v0) : std::make_tuple(alpha1, uei, v1);
            Q.emplace(bad_corner);

            return true;
    };

    auto project_vertex = [this] (const MatrixXd& V, const MatrixXi& /*F*/, const MatrixXi& /*E*/, 
            const MatrixXi& /*uE*/, const VectorXi& /*EMAP*/, const std::vector<std::vector<int>>& /*uE2E*/, 
            std::tuple<double,int,int> p, int e0, int e1, int /*v0*/, int /*v1*/) -> RowVector3d {
        
            // Projectee stored in tuple for easy retrival
            return utils::project_vertex_to_edge(V, e0, e1, std::get<2>(p));
    };

    auto pre_split = [this] (const MatrixXd& /*V*/, const MatrixXi& /*F*/, const MatrixXi& /*E*/, 
            const MatrixXi& /*uE*/, const VectorXi& /*EMAP*/, const std::vector<std::vector<int>>& /*uE2E*/, 
            std::priority_queue<std::tuple<double,int,int>> &/*Q*/, std::tuple<double,int,int> p, 
            int e0, int e1, int v0, int v1) -> bool {

            int projectee = std::get<2>(p);
            if (projectee != v0 && projectee != v1) {
                // Skip edge: triangle information is outdated (due to previous split)
                return false;
            }
            return true;
    };
    
    ArrayXb is_right_angled = ArrayXb::Constant(F_temp.rows(), false);
    auto post_split = [&bad_triangles, &is_right_angled, this] (const MatrixXd& V, const MatrixXi& F, const MatrixXi& E, 
            const MatrixXi& uE, const VectorXi& EMAP, const std::vector<std::vector<int>>& uE2E, 
            std::priority_queue<std::tuple<double,int,int>> &Q, int e0, int e1, int v0, 
            int v1, int vM, int f0, int f1, int f2, int f3) -> void {
        
            // Preallocation should sync with the internal preallocation in edge splitter
            if (V_poles.rows() != V.rows()) {
                assert(V_poles.rows() == vM);
                assert(is_split.size() == vM);
                assert(is_fixed.size() == vM);
                V_poles.conservativeResize(V.rows(), V_poles.cols());
                V_poles.bottomRows(V.rows() - vM).setConstant(-1);
                is_split.conservativeResize(V.rows());
                is_split.bottomRows(V.rows() - vM).setConstant(false);
                is_fixed.conservativeResize(V.rows());
                is_fixed.bottomRows(V.rows() - vM).setConstant(false);
            }

            // Determine projection side from orthogonality to split edge
            int projectee;
            // Check if {vM, v0} is ortho to {e0, e1}
            if (std::abs(utils::compute_internal_angle(V, e0, v0, vM) * 180 / M_PI - 90) < zero_TH) {
                assert(std::abs(utils::compute_internal_angle(V, e1, v0, vM) * 180 / M_PI - 90) < zero_TH);
                projectee = v0; // NOTE. f0, f2 are right angled
            }
            else { // Check if {vM, v1} is ortho to {e0, e1}
                assert (std::abs(utils::compute_internal_angle(V, e0, v1, vM) * 180 / M_PI - 90) < zero_TH);
                assert(std::abs(utils::compute_internal_angle(V, e1, v1, vM) * 180 / M_PI - 90) < zero_TH);
                projectee = v1; // NOTE. f1, f3 are right angled
            }

            // Project poles and mark split vertex
            V_poles.row(vM) = utils::project_vertex_to_edge(V_poles, e0, e1, projectee);
            is_split(vM) = true;

            // Check if non-projection split produced bad triangles and update queue
            std::array<int, 2> faces_to_check = (projectee == v0) ? std::array<int, 2>{f1, f3} : std::array<int, 2>{f0, f2};
            std::tuple<double, int, int>p;
            // for (int f : faces_to_check) {
            //     for (int c = 0; c < 3; c++) {
            //         int uei = EMAP(f + F.rows() * c);
            //         bool add_to_queue = bad_triangles(V, F, E, uE, EMAP, uE2E, Q, uei);
            //     }
            // }
            //  if (std::abs(utils::compute_internal_angle(V, e1, v1, vM) * 180 / M_PI - 90) > zero_TH){
            //         std::cout << "Error splitting face " << F.row(f2) << " at vertex " << vM << std::endl;
            //         std::cout << "Error splitting face " << F.row(f3) << " at vertex " << vM << std::endl;
            //         std::cout << "Input face: " << utils::compute_internal_angle(V, e0, e1, v0) * 180 / M_PI  << std::endl;
            //         std::cout << "Input face: " << utils::compute_internal_angle(V, e0, e1, v1) * 180 / M_PI  << std::endl;
            //         std::cout << is_split(e0) << " , " << is_split(e1) << " , " << is_split(v0) << " , " << is_split(v1) << std::endl;
            //     }
    };

    MatrixXd U;
    MatrixXi G;
    is_split.setConstant(false);
    bool clean_finish = utils::splitter(V_temp, F_temp, bad_triangles, project_vertex, pre_split, post_split, U, G);
    if (!clean_finish) {
        std::cerr << "Mesh not manifold after edge splitting." << std::endl;
    }

    // Audit topology changes
    int V_delta = U.rows() - V_temp.rows();
    int F_delta = G.rows() - F_temp.rows();
    assert(F_delta == 2 * V_delta);
    std::cout << "Increase in vertex count after splitting: " << V_delta << std::endl;
    std::cout << "Increase in face count after splitting: " << F_delta << std::endl;
    // std::cout << "Number of skipped splits" << skipped_count++ << std::endl;

    V_temp = U;
    F_temp = G;

    // Resize to remove any excess memory allocated
    int num_vertices = V_temp.rows();
    if (V_poles.rows() != num_vertices) {
        assert(V_poles.rows() > num_vertices);
        V_poles.conservativeResize(num_vertices, V_poles.cols());
        is_split.conservativeResize(num_vertices);
        is_fixed.conservativeResize(num_vertices);
    }
}