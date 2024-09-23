#pragma once

#include <Eigen/Core>
#include "utils.h"

class RemeshingTest;

class Remeshing {
    friend class RemeshingTest;

public:
    Remeshing(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, double scale = 0.95);
    
    void remesh();
    void setTargetEdgeLength(double scale);
   
    const Eigen::MatrixXd& getVertices() const { return V_remeshed; }
    const Eigen::MatrixXi& getFaces() const { return F_remeshed; }

private:
    /** \brief Split edges longer then 4/3 l at their midpoint. */
    void splitLongEdges();

    /** \brief Collapse edges shorter than 4/5 l at their midpoint. */
    void collapseShortEdges();

    /** \brief Flip edges to minimize the deviation from ideal valence (6) */
    void flipOptimalEdges();

    /** \brief Iteratively project gravity-weighted centroid to tangent plane */
    void tangentialSmoothing(double lambda = 1.0);

    const Eigen::MatrixXd& V; // Store reference to the original matrix
    const Eigen::MatrixXi& F;

    Eigen::MatrixXd V_remeshed;
    Eigen::MatrixXi F_remeshed;
    double target_edge_length, high_TH, low_TH; // target edge length (l)
};