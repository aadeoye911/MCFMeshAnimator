#pragma once

#include <vector>
#include <tuple>
#include <Eigen/Core>
#include <libqhullcpp/Qhull.h>

class MedialAxisTransform {
public:
    // Constructor: Initialize with the given points and faces
    MedialAxisTransform(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F);

    void compute();

    const Eigen::MatrixXd& getVertices() const { return V; }
    const Eigen::MatrixXi& getFaces() const { return F; }
    const Eigen::MatrixXd& getMedialPoles() const { return V_medial; }

private:

    bool isOutsideBoundingBox(Eigen::RowVector3d& point, Eigen::RowVector3d& AA, Eigen::RowVector3d& BB);
    
    Eigen::RowVector3d vectorizeQhullPoint(const orgQhull::QhullPoint& p);
    
    void computeVoronoiDiagram();
    
    void convertEigenToQhullPoints();

    std::tuple<double, double, double> vectorToTuple(const Eigen::RowVector3d& vec);

    const Eigen::MatrixXd& V; // Store reference to the original matrix
    const Eigen::MatrixXi& F;
    Eigen::MatrixXd N, V_medial;

    orgQhull::Qhull qhull_;       // Qhull object to compute Voronoi diagram
    std::vector<double> qhull_points_; // Stores the Qhull points

    struct TupleHash {
        template <typename T1, typename T2, typename T3>
        std::size_t operator() (const std::tuple<T1, T2, T3>& tuple) const {
            auto hash1 = std::hash<T1>{}(std::get<0>(tuple));
            auto hash2 = std::hash<T2>{}(std::get<1>(tuple));
            auto hash3 = std::hash<T3>{}(std::get<2>(tuple));
            return hash1 ^ (hash2 << 1) ^ (hash3 << 2); // Combine the hashes
        }
    };
};