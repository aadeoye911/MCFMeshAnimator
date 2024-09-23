
#include <Eigen/Core>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/unique_simplices.h>
#include <igl/readOBJ.h>
#include "remeshing.h"
#include <boost/filesystem.hpp>
#include "medial_axis_transform.h"
#include "mcf_skeleton.h"
#include "visualization.h"
#include "utils.h"

int main() {
    igl::opengl::glfw::Viewer viewer;
    
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

    std::string test_mesh = "../data/test_meshes/test_remeshed.obj";
    if (!igl::readOBJ(test_mesh, V, F)) {
        std::cerr << "Failed to load mesh from " << test_mesh << std::endl;
        return -1;  // Exit if the mesh cannot be loaded
    }

    // Check if the mesh is watertight
    if (!utils::is_watertight(F)) {
        std::cerr << "Input mesh must be watertight" << std::endl;
        return -1;  // Exit if the mesh is not watertight
    }
    
    // Remeshing remeshing(V, F);
    // remeshing.remesh();

    MedialAxisTransform mat(V, F);
    mat.compute();
    
    // Setup MCF Skeleton
    MCFSkeleton mcf(mat.getVertices(), mat.getFaces(), mat.getMedialPoles());

    int max_iterations = 20;
    mcf.skeletonize(5);
    view_skeleton(viewer, mcf);
    
    return 0;
}