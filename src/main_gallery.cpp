#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOBJ.h>
#include <igl/stb/write_image.h>
#include <boost/filesystem.hpp>
#include <Eigen/Core>
#include <iostream>

#include "remeshing.h"
#include "medial_axis_transform.h"
#include "mcf_skeleton.h"
#include "visualization.h"

namespace fs = boost::filesystem;

int main() {
    igl::opengl::glfw::Viewer viewer;

    std::string data_folder = "../data/meshes";
    std::string skel_folder = "../gallery/skeletons";
    std::string seg_folder = "../gallery/segmentations";

    // Ensure the output directories exist
    fs::create_directories(skel_folder);
    fs::create_directories(seg_folder);

    // Iterate over each OBJ file in the data folder
    for (const auto& entry : fs::directory_iterator(data_folder)) {
        if (entry.path().extension() == ".obj") {
            Eigen::MatrixXd V;
            Eigen::MatrixXi F;

            // Load the mesh
            if (!igl::readOBJ(entry.path().string(), V, F)) {
                std::cerr << "Failed to load mesh from " << entry.path() << std::endl;
                continue;
            }
            else{
                std::cout << "Loaded mesh: " << entry.path() << std::endl;
            }

            if (!utils::is_watertight(F)) {
                std::cerr << "Input mesh must be watertight" << std::endl;
                continue;  // Skip to next mesh
            }

            Remeshing remeshing(V, F);
            remeshing.remesh();
            
            MedialAxisTransform mat(remeshing.getVertices(), remeshing.getFaces());
            mat.compute();

            MCFSkeleton mcf(mat.getVertices(), mat.getFaces(), mat.getMedialPoles());
            mcf.skeletonize();

            // Allocate temporary buffers
            Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> R(1280,800);
            Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> G(1280,800);
            Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> B(1280,800);
            Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> A(1280,800);

            // Visualization 1: Skeleton
            view_skeleton(viewer, mcf);
            viewer.core().draw_buffer(viewer.data(),false,R,G,B,A);
            std::string skel_path = skel_folder + "/" + entry.path().stem().string() + "_skeleton.png";
            igl::stb::write_image(skel_path, R, G, B, A);

            // Visualization 2: Segmentation
            view_segmentation(viewer, mcf);
            viewer.core().draw_buffer(viewer.data(),false,R,G,B,A);
            std::string seg_path = seg_folder + "/" + entry.path().stem().string() + "_segmented.png";
            igl::stb::write_image(seg_path, R, G, B, A);
        }
    }
    return 0;
}