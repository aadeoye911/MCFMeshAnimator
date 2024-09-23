#pragma once

#include <Eigen/Core>
#include "mean_curvature_flow.h"
#include "utils.h"

class MCFSkeleton : public MeanCurvatureFlow {
public:
    MCFSkeleton(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const Eigen::MatrixXd& V_medial);

    void skeletonize(int max_interations = 50);  // Main skeletonization function
    
    void setParameters(double omega_L, double omega_H, double omega_M, double scale, double zero_TH, double alpha_TH);

    const Eigen::MatrixXd& V;
    const Eigen::MatrixXi& F;

    const Eigen::MatrixXd& getVertices() const { return V_temp; }
    const Eigen::MatrixXi& getFaces() const { return F_temp; }
    
    const Eigen::MatrixXd& getSkeletonVertices() const { return V_skeleton; }
    const Eigen::MatrixXi& getSkeletonFaces() const { return F_skeleton; }

    const Eigen::VectorXi& getCorrespondence() const {return corr; }
    int getIteration() const { return iter; }

private:
    void collapseShortEdges();
    void splitBadTriangles();
    void detectDegeneracies();

    Eigen::MatrixXd V_skeleton;
    Eigen::MatrixXi F_skeleton;
    Eigen::VectorXi corr;

    double edge_TH, alpha_TH; // remeshing constraints
    int iter;
};