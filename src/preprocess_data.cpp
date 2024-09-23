#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <boost/filesystem.hpp>
#include <Eigen/Core>
#include <iostream>
#include "remeshing.h"

namespace fs = boost::filesystem;

int main() {
    std::string data_folder = "../data/raw";
    std::string remeshed_folder = "../data/remeshed";

    // Ensure the output directories exist
    fs::create_directories(remeshed_folder);

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
            
            std::string output_path = remeshed_folder + "/" + entry.path().stem().string() + "_iso.obj";
            bool print = igl::writeOBJ(output_path, remeshing.getVertices(), remeshing.getFaces());
            if (!print) {
                std::cerr << "Failed to write output for " << entry.path() << std::endl;
            }
        }
    }
    return 0;
}