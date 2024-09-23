#include <iostream>
#include <limits>
#include <functional>
#include <cmath>
#include <unordered_set>
#include <Eigen/Core>

// libigl includes
#include <igl/decimate.h>
#include <igl/edge_lengths.h>
#include <igl/shortest_edge_and_midpoint.h>
#include <igl/decimate_callback_types.h>

#include "remeshing.h" // Class declaration and utility functions

using namespace Eigen;

void Remeshing::collapseShortEdges() {
    auto cost_and_placement = [this](int e, const MatrixXd &V, const MatrixXi &F, 
            const MatrixXi &E, const VectorXi &EMAP, const MatrixXi &EF, 
            const MatrixXi &EI, double &cost, RowVectorXd &p) -> void {
                
        igl::shortest_edge_and_midpoint(e, V, F, E, EMAP, EF, EI, cost, p);
        if (cost > low_TH){
            cost = std::numeric_limits<double>::infinity();
        }
    };

    // Stop decimation when the next edge is longer than threshold
    auto stopping_condition = [this](const MatrixXd &V, const MatrixXi &/*F*/,
            const MatrixXi &E, const VectorXi &/*EMAP*/, const MatrixXi &/*EF*/,
            const MatrixXi &/*EI*/, const igl::min_heap< std::tuple<double,int,int> > &Q,
            const VectorXi &/*EQ*/, const MatrixXd &/*C*/, const int /*e*/, const int /*e1*/,
            const int /*e2*/, const int /*f1*/, const int /*f2*/) -> bool {

        if (Q.empty()){
            return true;
        }
        auto next = Q.top(); // Access next edges
        int edge_index = std::get<1>(next); 
        double edge_length = (V.row(E(edge_index,0)) - V.row(E(edge_index, 1))).norm();
        
        return (edge_length > low_TH);
    };

    MatrixXd V_temp;
    MatrixXi F_temp;
    VectorXi J, I;
    bool clean_finish = igl::decimate(V_remeshed, F_remeshed, cost_and_placement, stopping_condition, V_temp, F_temp, J, I);
    if (!clean_finish) {
        std::cerr << "Decimation did not finish cleanly. Refine stopping condition." << std::endl;
    }

    int V_delta = V_remeshed.rows() - V_temp.rows();
    int F_delta = F_remeshed.rows() - F_temp.rows();
    assert(F_delta == 2 * V_delta);
    std::cout << "Decrease in vertex count after collapse: " << V_delta << std::endl;
    std::cout << "Decrease in face count after collapse: " << F_delta << std::endl;

    V_remeshed = V_temp;
    F_remeshed = F_temp;
}