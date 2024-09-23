#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>

using ArrayXb = Eigen::Array<bool, Eigen::Dynamic,1>;

class MeanCurvatureFlow {
public:
    MeanCurvatureFlow(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const Eigen::MatrixXd& V_medial);

protected:
    void updateLaplacian();
    void updateWeights();
    void performMeshContraction();
    virtual void setParameters(double omega_L, double omega_H, double omega_M, double zero_TH);
    
    const Eigen::MatrixXd& V, V_medial; // Store reference to the original matrix
    const Eigen::MatrixXi& F;

    Eigen::MatrixXd V_temp, V_poles;
    Eigen::MatrixXi F_temp;

    Eigen::SparseMatrix<double> L; // Laplacian and weight matrices
    Eigen::VectorXd W_L, W_H, W_M;
    double omega_L, omega_H, omega_M, zero_TH; // Constraint parameters
    ArrayXb is_fixed, is_split;
};