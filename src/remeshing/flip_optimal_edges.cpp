#include <iostream>
#include <vector>
#include <cassert>
#include <algorithm>
#include <utility> 
#include <Eigen/Core>

// libigl includes
#include <igl/unique_edge_map.h>
#include <igl/adjacency_list.h>
#include <igl/flip_edge.h>

#include "remeshing.h" // Class declaration and utility functions

using namespace Eigen;

using ArrayXb = Eigen::Array<bool, Eigen::Dynamic, 1>;

void Remeshing::flipOptimalEdges() {
    int num_faces = F_remeshed.rows();
    int num_vertices = F_remeshed.maxCoeff() + 1;
    MatrixXi E, uE;
    VectorXi EMAP;
    std::vector<std::vector<int>> uE2E, A;
    igl::unique_edge_map(F_remeshed, E, uE, EMAP, uE2E);

    // Compute vertex valences for mesh
    VectorXi current_vals(num_vertices);
    igl::adjacency_list(F_remeshed, A);
    for (int v = 0; v < num_vertices; v++){
        current_vals[v] = A[v].size();
    }
    double initial_deviation = (current_vals.array() - 6).cwiseAbs().sum();

    VectorXi N_edge(4), delta_vals(4), target_vals(4);
    delta_vals << -1, -1, 1, 1; // Valence changes: -1 for edge vertices and +1 for opposite vertices
    target_vals << 6, 6, 6, 6; // Optimal valence of 6 assuming watertight mesh

    ArrayXb is_flipped = ArrayXb::Constant(uE.rows(), false);       
    for (int uei = 0; uei < uE.rows(); uei++) {

        if(uE2E[uei].size() != 2) {continue;}

        const size_t f1 = uE2E[uei][0] % num_faces;
        const size_t f2 = uE2E[uei][1] % num_faces;
        const size_t c1 = uE2E[uei][0] / num_faces;
        const size_t c2 = uE2E[uei][1] / num_faces;
        assert(c1 < 3);
        assert(c2 < 3);
        assert(f1 != f2);

        // Define edge neighbourhood as endpoints and opposite vertices
        const size_t v1 = F_remeshed(f1, (c1+1)%3);
        const size_t v2 = F_remeshed(f1, (c1+2)%3);
        const size_t v4 = F_remeshed(f1, c1);
        const size_t v3 = F_remeshed(f2, c2);
        assert(F_remeshed(f2, (c2+2)%3) == v1);
        assert(F_remeshed(f2, (c2+1)%3) == v2);
    
        N_edge << v1, v2, v3, v4;

        // Compute deviation before and after flipping
        double deviation_before = (current_vals(N_edge) - target_vals).cwiseAbs().sum();
        double deviation_after = (current_vals(N_edge) + delta_vals - target_vals).cwiseAbs().sum();

        if (deviation_after < deviation_before) {
            // Skip non-manifold flips
            if(std::count(A[v3].begin(),A[v3].end(),v4)){continue;}

            igl::flip_edge(F_remeshed, E, uE, EMAP, uE2E, uei);
            current_vals(N_edge) += delta_vals;
            is_flipped[uei] = true;

            // Update adjacency list
            A[v1].erase(std::remove_if(A[v1].begin(), A[v1].end(), [&v2](const int & v) { return v == v2; }), A[v1].end());
            A[v2].erase(std::remove_if(A[v2].begin(), A[v2].end(), [&v1](const int & v) { return v == v1; }), A[v2].end());
            A[v3].push_back(v4);
            A[v4].push_back(v3);
        }
    }

    double final_deviation = (current_vals.array() - 6).cwiseAbs().sum();
    std::cout << "Flipped count: " << is_flipped.count() << ". Total edges: " << uE.rows() << std::endl;
    std::cout << "Total deviation before: " << initial_deviation << ". Total deviation after: " << final_deviation << std::endl;

    assert(utils::is_watertight(F_remeshed));
}