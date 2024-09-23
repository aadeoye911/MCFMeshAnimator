#include <Eigen/Core>
#include <igl/per_vertex_normals.h>
#include "visualization.h"

void view_remeshing(igl::opengl::glfw::Viewer &viewer, const Remeshing &remesh) {
    viewer.data().clear();  // Clear the viewer to prepare for new visualization

    // Get the vertices and faces from the remeshing class
    const Eigen::MatrixXd &V = remesh.getVertices();
    const Eigen::MatrixXi &F = remesh.getFaces();

    // Set the remeshed mesh
    viewer.data().set_mesh(V, F);
    viewer.data().show_lines = true;
    viewer.data().line_width = 1.5f;
    viewer.launch();
}

void view_mat(igl::opengl::glfw::Viewer &viewer, const MedialAxisTransform &mat) {
    viewer.data().clear();

    const Eigen::MatrixXd &V = mat.getVertices();
    const Eigen::MatrixXi &F = mat.getFaces();

    viewer.data().set_mesh(V, F);
    viewer.data().set_colors(Eigen::RowVector3d(1, 1, 1)); // White color
    viewer.data().show_faces = false; // Make faces transparent
    viewer.data().show_overlay_depth = false;

    viewer.data().add_points(mat.getMedialPoles(), Eigen::RowVector3d(0, 1, 0)); // Green points for MAT
    viewer.launch();
}

void view_skeleton(igl::opengl::glfw::Viewer &viewer, const MCFSkeleton &mcf) {
    viewer.data().clear();
    viewer.data().set_mesh(mcf.V, mcf.F);
    viewer.data().set_face_based(true);
    viewer.data().show_faces = false;
    viewer.data().show_lines = true;
    
    // Set colors with transparency (RGBA)
    Eigen::MatrixXd C(mcf.F.rows(), 4); // Color matrix with alpha channel
    C.col(0) = Eigen::VectorXd::Constant(mcf.F.rows(), 1.0); // Red channel
    C.col(1) = Eigen::VectorXd::Constant(mcf.F.rows(), 1.0); // Green channel
    C.col(2) = Eigen::VectorXd::Constant(mcf.F.rows(), 1.0); // Blue channel
    C.col(3) = Eigen::VectorXd::Constant(mcf.F.rows(), 0.2); // Alpha channel (0.5 for 50% transparency)
    
    viewer.data().set_colors(C);
    viewer.core().background_color.setOnes(); // Set background to white for visibility

    int skeleton_view = viewer.append_mesh();
    viewer.data(skeleton_view).set_colors(Eigen::RowVector3d(1, 0, 0)); // Red color for MCF skeleton
    viewer.data(skeleton_view).set_face_based(true);
    viewer.data(skeleton_view).show_faces = true;
    viewer.data(skeleton_view).set_mesh(mcf.getVertices(), mcf.getFaces());

    viewer.launch(); 
}

void view_segmentation(igl::opengl::glfw::Viewer &viewer, const MCFSkeleton &mcf) {
    viewer.data().clear();
    viewer.data().set_mesh(mcf.V, mcf.F);
    Eigen::MatrixXd N;
    igl::per_vertex_normals(mcf.getSkeletonVertices(), mcf.getSkeletonFaces(), N);

    Eigen::MatrixXd C;
    C = N(mcf.getCorrespondence(), Eigen::all).rowwise().normalized();
    viewer.data().set_colors(C);
    viewer.core().background_color.setOnes(); // Set background to white for visibility
    viewer.launch();
}

