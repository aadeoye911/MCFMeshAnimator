#include <iostream>
#include <vector>
#include <limits>
#include <Eigen/Core>

// libigl includes
#include <igl/adjacency_list.h>
#include <igl/decimate.h>
#include <igl/decimate_trivial_callbacks.h>
#include <igl/edge_lengths.h>
#include <igl/infinite_cost_stopping_condition.h>
#include <igl/shortest_edge_and_midpoint.h>

#include "mcf_skeleton.h" // Subclass declaration and utility functions

using namespace Eigen;
using ArrayXb = Eigen::Array<bool, Eigen::Dynamic,1>;

void MCFSkeleton::collapseShortEdges() {
    // Collapse shortest edges at their midpoint
    auto cost_and_placement = [this](int e, const MatrixXd &V, const MatrixXi &F, const MatrixXi &E, const VectorXi &EMAP, 
                            const MatrixXi &EF, const MatrixXi &EI, double &cost, RowVectorXd &p) -> void {
        // Default cost and placement strategy assigns cost as equal to edge length
        igl::shortest_edge_and_midpoint(e, V, F, E, EMAP, EF, EI, cost, p);
        // Reassign infinite cost to edge longer than threshold
        if (cost > edge_TH) {
            cost = std::numeric_limits<double>::infinity();
        }
        // Do not collapse fixed edges
        if (is_fixed(E(e, 0)) && is_fixed(E(e, 1))) {
            cost = std::numeric_limits<double>::infinity();
        }
    };

    // Pre collapse: trivial callback to always attempt collapse
    int s, d; // Store original edge index
    auto pre_collapse = [&](const MatrixXd &V, const MatrixXi &F, const MatrixXi &E, const VectorXi &EMAP, 
                        const MatrixXi &EF, const MatrixXi &EI, const igl::min_heap<std::tuple<double,int,int>> &/*Q*/,
                        const VectorXi &/*EQ*/, const MatrixXd &/*C*/, const int e) -> bool {
        // Log degenerate edge vertices for post collapse
        const int eflip = E(e, 0) > E(e, 1);
        s = eflip ? E(e, 1) : E(e, 0);  // Index of retained vertex
        d = eflip ? E(e, 0) : E(e, 1);  // Index of collapsed vertex

        return true;  // Allow the collapse
    };

    // Post collapse: Update correspondence and poles
    ArrayXb is_null_vertex = ArrayXb::Constant(V_temp.rows(), false);
    auto post_collapse = [&](const MatrixXd &V, const MatrixXi &F, const MatrixXi &/*E*/, const VectorXi &/*EMAP*/,
                        const MatrixXi &/*EF*/, const MatrixXi &/*EI*/, const igl::min_heap<std::tuple<double,int,int>> &/*Q*/,
                        const VectorXi &/*EQ*/, const MatrixXd &/*C*/, const int /*e*/, const int /*e1*/, const int /*e2*/,
                        const int /*f1*/, const int /*f2*/, const bool collapsed) -> void {
        if (collapsed) {
            assert (V.row(s) == V.row(d));
            is_null_vertex(d) = true; // Mark collapsed vertex 
            corr = (corr.array() == d).select(s, corr); // Update correspondence

            // Update closest pole for retained vertex
            if ((V_poles.row(s) - V.row(s)).squaredNorm() > (V_poles.row(d) - V.row(s)).squaredNorm()) {
                V_poles.row(s) = V_poles.row(d);
            }
        }
    };

    // Stop decimation if all remeaining edges have infinite cost
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
        
        return (edge_length > edge_TH);
    };

    MatrixXd U;
    MatrixXi G;
    VectorXi J, I;
    bool clean_finish = igl::decimate(V_temp, F_temp, cost_and_placement, stopping_condition, pre_collapse, post_collapse, U, G, J, I);
    if (!clean_finish) {
        std::cerr << "Decimation did not finish cleanly. Refine stopping condition" << std::endl;
    }
    assert(utils::is_watertight(F_temp));

    int V_delta = V_temp.rows() - U.rows();
    int F_delta = F_temp.rows() - G.rows();
    assert(F_delta == 2 * V_delta);
    std::cout << "Decrease in vertex count after collapse: " << V_delta << std::endl;
    std::cout << "Decrease in face count after collapse: " << F_delta << std::endl;
    std::cout << "Number of vertices after collapse: " << V_temp.rows() << std::endl;

    V_temp = U;
    F_temp = G;

    ArrayXb is_fixed_filtered(V_temp.rows());
    ArrayXb is_split_filtered(V_temp.rows());
    MatrixXd V_poles_filtered(V_temp.rows(), V_temp.cols());

    int e = 0;
    for (int i = 0; i < is_null_vertex.size(); ++i) {
        if (!is_null_vertex(i)) {
            V_poles_filtered.row(e) = V_poles.row(i);
            is_fixed_filtered(e) = is_fixed(i);
            is_split_filtered(e) = is_split(i);
            e++;
        }
    }
    V_poles = V_poles_filtered;
    is_split = is_split_filtered;
    is_fixed = is_fixed_filtered;
}