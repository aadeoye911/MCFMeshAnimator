#include <limits>
#include <unordered_map>
#include <Eigen/Core>

// libigl includes
#include <igl/per_vertex_normals.h>

// Qhull includes
#include <libqhullcpp/QhullQh.h>
#include <libqhullcpp/QhullPoint.h>
#include <libqhullcpp/QhullFacet.h>
#include <libqhullcpp/QhullFacetSet.h>
#include <libqhullcpp/QhullLinkedList.h>
#include <libqhullcpp/QhullVertex.h>
#include <libqhullcpp/QhullVertexSet.h>

#include "medial_axis_transform.h" // Class declaration

using namespace Eigen;
using namespace orgQhull;

using ArrayXb = Eigen::Array<bool, Eigen::Dynamic,1>;

// Constructor: Initialize with the given points
MedialAxisTransform::MedialAxisTransform(const MatrixXd& V, const MatrixXi& F)
    : V(V), F(F) {
    igl::per_vertex_normals(V, F, N);
    convertEigenToQhullPoints();
    computeVoronoiDiagram(); 
}

void MedialAxisTransform::compute() {
    V_medial.resize(V.rows(), 3);
    RowVector3d min_coord = V.colwise().minCoeff(); // lower bounding coord
    RowVector3d max_coord = V.colwise().maxCoeff(); // upper bounding coord
    RowVector3d furthest_pole, voronoi_vertex;

    // Map indexing Qhull points to V matrix indices
    std::unordered_map<std::tuple<double, double, double>, int, TupleHash> coord_to_index_map;
    for (int i = 0; i < V.rows(); ++i) {
        coord_to_index_map[vectorToTuple(V.row(i))] = i;
    }

    QhullVertexListIterator v(qhull_.vertexList());
    while(v.hasNext()) {
        QhullVertex input_vertex = v.next();
        RowVector3d input_vec = vectorizeQhullPoint(input_vertex.point());
        auto vertex_finder = vectorToTuple(input_vec);
        
        if (coord_to_index_map.find(vertex_finder) != coord_to_index_map.end()) {
            int v_index = coord_to_index_map[vertex_finder];

            // Iterate over facets to find furthest Voronoi pole
            double dist_prev = std::numeric_limits<double>::infinity();
            bool found_pole = false; // flag to track whether valid pole was found
    
            QhullFacetSetIterator f(input_vertex.neighborFacets());
            while(f.hasNext()) {
                QhullFacet facet = f.next();
                if (facet.isUpperDelaunay()) { continue; } // Ignore upper facets of convex hull

                voronoi_vertex = vectorizeQhullPoint(facet.voronoiVertex());
                if (isOutsideBoundingBox(voronoi_vertex, min_coord, max_coord)) { continue; } // Ignore poles outside bounding box

                double dist = N.row(v_index).dot((voronoi_vertex - V.row(v_index)));
                if (dist < 0 && dist < dist_prev) {
                    dist_prev = dist;
                    furthest_pole = voronoi_vertex;
                    found_pole = true;
                }
            }
            if (found_pole) {
                V_medial.row(v_index) = furthest_pole;
            }
        }
    }
    assert(V_medial.rowwise().norm().allFinite());
}

bool MedialAxisTransform::isOutsideBoundingBox(RowVector3d& point, RowVector3d& AA, RowVector3d& BB) {
    return !((point.array() >= AA.array()).all() && (point.array() <= BB.array()).all());
}

RowVector3d MedialAxisTransform::vectorizeQhullPoint(const QhullPoint& p) {
    return RowVector3d(p[0], p[1], p[2]);
}

std::tuple<double, double, double> MedialAxisTransform::vectorToTuple(const RowVector3d& vec) {
    return std::make_tuple(vec[0], vec[1], vec[2]);
}

void MedialAxisTransform::computeVoronoiDiagram(){
    qhull_.runQhull("", 3, qhull_points_.size() / 3, qhull_points_.data(), "v Qbb");
}

void MedialAxisTransform::convertEigenToQhullPoints() {
    // Tranpose to make vertex matrix column major
    Eigen::MatrixXd V_transposed = V.transpose();
    qhull_points_.assign(V_transposed.data(), V_transposed.data() + V_transposed.size());
}