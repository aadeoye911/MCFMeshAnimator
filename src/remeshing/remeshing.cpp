#include <iostream>
#include <igl/avg_edge_length.h>
#include "remeshing.h" // Class declaration and utility functions

using namespace Eigen;

Remeshing::Remeshing(const MatrixXd& V, const MatrixXi& F, double scale) 
    : V(V), F(F), V_remeshed(V), F_remeshed(F) {
    setTargetEdgeLength(scale);
}

void Remeshing::remesh() {
    std::cout << "Vertex count before remeshing: " << V.rows() << std::endl; 
    std::cout << "Face count before remeshing: " << F.rows() << std::endl; 

    for (int iter = 0; iter < 20; ++iter) {
        splitLongEdges();
        collapseShortEdges();
        flipOptimalEdges();
        tangentialSmoothing();
    }
    std::cout << "Vertex count after remeshing: " << V_remeshed.rows() << std::endl; 
    std::cout << "Face count after remeshing: " << F_remeshed.rows() << std::endl; 
}

void Remeshing::setTargetEdgeLength(double scale) {
    target_edge_length = scale * igl::avg_edge_length(V, F);
    high_TH = 4.0 / 3.0 * target_edge_length;
    low_TH = 4.0 / 5.0 * target_edge_length;
}