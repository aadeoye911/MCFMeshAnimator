#include <gtest/gtest.h>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/opengl/glfw/Viewer.h>

#include "remeshing.h"
#include "medial_axis_transform.h"
#include "visualizaton.h"
#include "utils.h"

class MedialAxisTransformTest : public ::testing::Test {
protected:
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    
    virtual void SetUp() {
        std::string test_mesh = "data/test_meshes/test_mesh.obj";
        if (!igl::readOBJ(test_mesh, V, F)) {
            FAIL() << "Failed to load mesh from " << test_mesh;
        }

        // Check mesh is watertight
        ASSERT_EQ(utils::is_watertight(F), true) << "Test mesh is not watertight";
    }
};

TEST_F(MedialAxisTransformTest, MedialAxisWithoutRemeshing) {
    MedialAxisTransform mat(V, F);
    mat.compute();

    igl::viewer::Viewer viewer;
    view_mat(viewer, mat);
}

TEST_F(MedialAxisTransformTest, MedialAxisWithRemeshing) {
    Remeshing remesh(V, F);
    MedialAxisTransform mat(remesh.getVertices(), remesh.getFaces());
    mat.compute();

    igl::opengl::glfw::Viewer viewer;
    view_mat(viewer, mat);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}