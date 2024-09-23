#include <functional>
#include <cmath>
#include <Eigen/Core>

// libigl includes
#include <igl/is_edge_manifold.h>
#include <igl/edge_flaps.h>

#include "remeshing.h" // Class declaration and utility functions

using namespace Eigen;

void Remeshing::splitLongEdges() {    
    auto long_edges = [this] (const MatrixXd& V, const MatrixXi& /*F*/, const MatrixXi& uE, 
            const MatrixXi& /*uE*/, const VectorXi& /*EMAP*/, const std::vector<std::vector<int>>& /*uE2E*/, 
            std::priority_queue<std::tuple<double,int,int>> &Q, int uei) -> bool {
        
        int e0 = uE(uei, 0);
        int e1 = uE(uei, 1);
        double length = (V.row(e0) - V.row(e1)).norm();
        if (length > high_TH) {
            Q.emplace(length, uei, 0);
            return true;
        }
        return false;
    };

    auto midpoint_vertex = [] (const MatrixXd& V, const MatrixXi& /*F*/, const MatrixXi& /*E*/, 
            const MatrixXi& /*uE*/, const VectorXi& /*EMAP*/, const std::vector<std::vector<int>>& /*uE2E*/, 
            std::tuple<double,int,int>&/*p*/, int e0, int e1, int /*v0*/, int /*v1*/) -> RowVector3d {
            
        // Assign split vertex at edge midpoint
        return (V.row(e0) + V.row(e1)) / 2.0;
    };

    // Trivial callback: edge length is a constant condition
    auto pre_split = [] (const MatrixXd& /*V*/, const MatrixXi& /*F*/, const MatrixXi& /*E*/, 
            const MatrixXi& /*uE*/, const VectorXi& /*EMAP*/, const std::vector<std::vector<int>>& /*uE2E*/, 
            std::priority_queue<std::tuple<double,int,int>>& /*Q*/, std::tuple<double,int,int>& /*p*/, 
            int /*e0*/, int /*e1*/, int /*v0*/, int /*v1*/) -> bool {
            
        return true;
    };

    // Trivial callback: long split edges are not added to the queue
    auto post_split = [] (const MatrixXd& /*V*/, const MatrixXi& /*F*/, const MatrixXi& /*E*/, 
            const MatrixXi& /*uE*/, const VectorXi& /*EMAP*/, const std::vector<std::vector<int>>& /*uE2E*/, 
            std::priority_queue<std::tuple<double,int,int>> &/*Q*/, int /*e0*/, int /*e1*/, int /*v0*/, 
            int /*v1*/, int /*vM*/, int /*f0*/, int /*f1*/, int /*f2*/, int /*f3*/) -> void {

        return;
    };

    MatrixXd V_temp;
    MatrixXi F_temp;
    bool clean_finish = utils::splitter(V_remeshed, F_remeshed, long_edges, midpoint_vertex, pre_split, post_split, V_temp, F_temp);
    if (!clean_finish) {
        std::cerr << "Mesh not manifold after edge splitting." << std::endl;
    }

    int V_delta = V_remeshed.rows() - V_temp.rows();
    int F_delta = F_remeshed.rows() - F_temp.rows();
    assert(F_delta == 2 * V_delta);
    std::cout << "Increase in vertex count after collapse: " << V_delta << std::endl;
    std::cout << "Increase in face count after collapse: " << F_delta << std::endl;

    V_remeshed = V_temp;
    F_remeshed = F_temp;
}