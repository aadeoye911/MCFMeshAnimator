#include <gtest/gtest.h>
#include <Eigen/Core>

// libigl includes
#include <igl/readOBJ.h>
#include <igl/edge_lengths.h>
#include <igl/opengl/glfw/Viewer.h>

#include "remeshing.h"
#include "utils.h"
#include "visualization.h"

class RemeshingTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Load the mesh from file into V and F
        std::string test_mesh = "data/test_meshes/test_mesh.obj";
        if (!igl::readOBJ(test_mesh, V, F)) {
            FAIL() << "Failed to load mesh from " << test_mesh;
        }

        // Check mesh is watertight
        ASSERT_EQ(utils::is_watertight(F), true) << "Test mesh is not watertight";
        
        // Initialize the Remeshing object
        remesh = new Remeshing(V, F);
    }

    void TearDown() override {
        // Clean up the Remeshing object
        delete remesh;
    }
    isotropicRemeshing* remesh;
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
};

TEST_F(RemeshingTest, checkFaceAndVertexCountAfterSplit) {
    std::cout << "Initial vertex count (before splitting): " << V.rows() << std::endl;
    std::cout << "Initial face count (before splitting): " << F.rows() << std::endl;

    double max_length = 4.0 / 3.0 * remesh->target_edge_length;
    Eigen::MatrixXd L;
    igl::edge_lengths(V, F, L);
    int edges_to_split = (L.array() > max_length).count();
    std::cout << "Number of edges identified for splitting: " << edges_to_split << std::endl;
    
    remesh->split_long_edges();

    // Assert that the face count increase equals twice the increase in vertex count
    ASSERT_EQ(remesh->F_temp.rows()- F.rows(), 2 * (remesh->V_temp.rows() - V.rows()))
        << "Face count increase should equal twice the increase in vertex count.";
}

TEST_F(RemeshingTest, checkEdgeLengthsAfterSplit) {
    double max_length = 4.0 / 3.0 * remesh->target_edge_length;
    std::cout << "Max length threshold: " << max_length << std::endl;

    Eigen::MatrixXd L;
    igl::edge_lengths(V, F, L);
    double longest_edge_before = L.maxCoeff();
    std::cout << "Initial longest edge length: " << longest_edge_before << std::endl;
    
    remesh->split_long_edges();

    igl::edge_lengths(remesh->V_temp, remesh->F_temp, L);
    double longest_edge_after = L.maxCoeff();
    std::cout << "Final longest edge length: " << longest_edge_after << std::endl;

    ASSERT_LE(longest_edge_after, longest_edge_before)
        << "Longest edge length should decrease or remain the same after splitting.";
        
    ASSERT_LE(longest_edge_after, max_length)
        << "Longest edge length (after splitting) should be less than or equal to max length threshold.";    
}

TEST_F(RemeshingTest, checkFaceAndVertexCountAfterCollapse) {
    std::cout << "Initial vertex count (before collapse): " << V.rows() << std::endl;
    std::cout << "Initial face count (before collapse): " << F.rows() << std::endl;

    double min_length = 4.0 / 5.0 * remesh->target_edge_length;
    Eigen::MatrixXd L;
    igl::edge_lengths(V, F, L);
    int edges_to_collapse = (L.array() > min_length).count();
    std::cout << "Number of edges identified for collapse: " << edges_to_collapse << std::endl;
    
    remesh->collapse_short_edges();

    std::cout << "Final vertex count (after splitting): " << remesh->V_temp.rows()<< std::endl;
    std::cout << "Final face count (after splitting): " << remesh->F_temp.rows() << std::endl;

    // Assert that the face count increase equals twice the increase in vertex count
    ASSERT_EQ(F.rows() - remesh->F_temp.rows(), 2 * (V.rows() - remesh->V_temp.rows()))
        << "Face count decrease should equal twice the decrease in vertex count.";
    
    // Assert that the number of edges identified for splitting equals the number of new vertices added
    ASSERT_EQ(edges_to_collapse, V.rows() - remesh->V_temp.rows())
        << "Number of edges to collapse should equal number of degenerate vertices added.";
}

TEST_F(RemeshingTest, checkEdgeLengthsAfterCollapse) {
    double min_length = 4.0 / 5.0 * remesh->target_edge_length;
    std::cout << "Min length threshold: " << min_length << std::endl;

    Eigen::MatrixXd L;
    igl::edge_lengths(V, F, L);
    double shortest_edge_before = L.minCoeff();
    std::cout << "Initial shortest edge length: " << shortest_edge_before << std::endl;
    
    remesh->collapse_short_edges();

    igl::edge_lengths(remesh->V_temp, remesh->F_temp, L);
    double shortest_edge_after = L.minCoeff();
    std::cout << "Final shortest edge length: " << shortest_edge_after << std::endl;

    ASSERT_GE(shortest_edge_after, shortest_edge_before)
        << "Shortest edge length should increase or remain the same after collapse.";

    ASSERT_GE(shortest_edge_after, min_length)
        << "Shortest edge length (after splitting) should be greater than or equal to min length threshold.";    
}

TEST_F(RemeshingTest, testFullRemesh) {
    igl::opengl::glfw::Viewer viewer;
    remesh.remesh();
    view_remeshing(viewer, remesh);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}