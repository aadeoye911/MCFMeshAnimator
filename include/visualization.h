#include <igl/opengl/glfw/Viewer.h>
#include "remeshing.h"
#include "medial_axis_transform.h"
#include "mcf_skeleton.h"

void view_remeshing(igl::opengl::glfw::Viewer &viewer, const Remeshing &remesh);

void view_mat(igl::opengl::glfw::Viewer &viewer, const MedialAxisTransform &mat);

void view_skeleton(igl::opengl::glfw::Viewer &viewer, const MCFSkeleton &mcf);

void view_segmentation(igl::opengl::glfw::Viewer &viewer, const MCFSkeleton &mcf);

