#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Sparse>

// libigl includes
#include <igl/cat.h>
#include <igl/cotmatrix.h>
#include <igl/invert_diag.h>
#include <igl/massmatrix.h>
#include <igl/adjacency_matrix.h>

#include "mean_curvature_flow.h" // Base class 

using namespace Eigen;

using ArrayXb = Eigen::Array<bool, Eigen::Dynamic,1>;

MeanCurvatureFlow::MeanCurvatureFlow(const MatrixXd& V, const MatrixXi& F, const MatrixXd& V_medial) 
    : V(V), F(F), V_medial(V_medial), V_temp(V), F_temp(F), V_poles(V_medial) { 

    is_fixed = ArrayXb::Constant(V.rows(), false); // tracking variable
    is_split = ArrayXb::Constant(V.rows(), false); // tracking variable
}

void MeanCurvatureFlow::performMeshContraction(){
    SparseMatrix<double> A (3*L.rows(), L.cols());
    A.reserve(L.nonZeros() + 2*L.rows());
    std::cout << "Start" << std::endl;
    for (int v = 0; v < L.cols(); ++v) {
        for (SparseMatrix<double>::InnerIterator it(L, v); it; ++it) {
            int row = it.row();
            int col = it.col();
            A.coeffRef(row, col) = W_L(row) * L.coeffRef(row, col);
            if (row == col) {
                A.coeffRef(row + L.rows(), col) = W_H(row); // Velocity constraint
                A.coeffRef(row + 2 * L.rows(), col) = W_M(row); // Pole constraint
            }
        }
    }
    A.makeCompressed();

    MatrixXd b(3 * L.rows(), 3); 
    b.middleRows(L.rows(), L.rows()) = W_H.asDiagonal() * V_temp;
    b.bottomRows(L.rows()) = W_M.asDiagonal() * V_poles;
    
    SparseMatrix<double> At = A.transpose();
    SimplicialLDLT<SparseMatrix<double>> solver;
    solver.compute(At * A);
    if (solver.info() != Eigen::Success) {
        std::cerr << "Solver failed!" << std::endl;
        return;
    }

    V_temp = solver.solve(At * b);
    assert(V_temp.allFinite());
}

void MeanCurvatureFlow::updateLaplacian() {
    SparseMatrix<double> area, adj, weights, C;
    igl::cotmatrix(V_temp, F_temp, C);
    igl::massmatrix(V_temp, F_temp, igl::MASSMATRIX_TYPE_VORONOI, area);
    igl::adjacency_matrix(F_temp, adj);

    weights = area * adj;
    VectorXd norm = 2 * weights.transpose() * VectorXd::Ones(weights.cols());
    L = norm.cwiseInverse().asDiagonal() * C;
}

void MeanCurvatureFlow::updateWeights(){
    W_L.resize(V_temp.rows()), W_H.resize(V_temp.rows()), W_M.resize(V_temp.rows());
    for (int i = 0; i < V_temp.rows(); ++i) {
        W_L(i) = is_fixed(i) ? 0.0 : omega_L;
        W_H(i) = is_fixed(i) ? 1.0/zero_TH : omega_H;
        W_M(i) = (is_fixed(i) || is_split(i)) ? 0.0 : omega_M;
    }    
}

void MeanCurvatureFlow::setParameters(double omega_L, double omega_H, double omega_M, double zero_TH) {
    this->omega_L = omega_L;
    this->omega_H = omega_H;
    this->omega_M = omega_M;
    this->zero_TH = zero_TH;
}