#include <iostream>
#include <Eigen/Core>
#include <Eigen/Sparse>

// libigl includes
#include <igl/adjacency_matrix.h>
#include <igl/massmatrix.h>
#include <igl/per_vertex_normals.h>
#include <igl/dot_row.h>

#include "remeshing.h"

using namespace Eigen;

void Remeshing::tangentialSmoothing(double lambda) {
    SparseMatrix<double> area, adj, weights;
    MatrixXd gravity, D, N;
    igl::adjacency_matrix(F_remeshed, adj);
    igl::massmatrix(V_remeshed, F_remeshed, igl::MASSMATRIX_TYPE_VORONOI, area);
    weights = area * adj;

    VectorXd norm(weights.rows());
    for (int k = 0; k < weights.outerSize(); ++k) {
        double colwise_sum = 0.0;
        for (SparseMatrix<double>::InnerIterator it(weights, k); it; ++it) {
            colwise_sum += it.value();
        }
        norm[k] = colwise_sum;
    }

    // Area-weighted gravity centroid
    gravity = norm.cwiseInverse().asDiagonal() * weights.transpose() * V_remeshed;

    // Project to tangential plane
    D = gravity - V_remeshed;
    igl::per_vertex_normals(V_remeshed, F_remeshed, N);
    VectorXd dot_product = igl::dot_row(N, D);
    V_remeshed += lambda * (D - MatrixXd(dot_product.replicate(1, 3).array() * N.array()));
}