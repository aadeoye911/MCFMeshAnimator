#pragma once

#include <iostream>
#include <stdexcept>
#include <vector>
#include <queue>
#include <algorithm>
#include <functional>
#include <Eigen/Core>
#include <igl/boundary_loop.h>
#include <igl/dot_row.h>
#include <igl/is_edge_manifold.h>
#include <igl/is_vertex_manifold.h>
#include <igl/unique_edge_map.h>
#include <igl/unique_rows.h>

namespace utils {

    /** \brief Check if the mesh is watertight (i.e., no holes or gaps). */
    using split_condition_callback = std::function<bool(const Eigen::MatrixXd& /*V*/, const Eigen::MatrixXi& /*F*/, 
                    const Eigen::MatrixXi &/*E*/, const Eigen::MatrixXi &/*uE*/, const Eigen::VectorXi& /*EMAP*/,
                    const std::vector<std::vector<int>> &/*uE2E*/, std::priority_queue<std::tuple<double,int,int>> &/*Q*/, int/*uei*/)>;

    /** \brief Check if the mesh is watertight (i.e., no holes or gaps). */
    using split_placement_callback = std::function<Eigen::RowVector3d(const Eigen::MatrixXd& /*V*/, const Eigen::MatrixXi&/*F*/, 
                    const Eigen::MatrixXi &/*E*/, const Eigen::MatrixXi &/*uE*/, const Eigen::VectorXi& /*EMAP*/,
                    const std::vector<std::vector<int>> &/*uE2E*/, std::tuple<double,int,int>& /*p*/, 
                    int/*e0*/, int/*e1*/, int /*v0*/, int /*v1*/)>;

    /** \brief Check if the mesh is watertight (i.e., no holes or gaps). */
    using pre_split_callback = std::function<bool(const Eigen::MatrixXd &/*V*/, const Eigen::MatrixXi &/*F*/, 
                    const Eigen::MatrixXi &/*E*/, const Eigen::MatrixXi &/*uE*/, const Eigen::VectorXi& /*EMAP*/,
                    const std::vector<std::vector<int>> &/*uE2E*/, std::priority_queue<std::tuple<double,int,int>> &/*Q*/, 
                    std::tuple<double,int,int>& /*p*/, int/*e0*/, int/*e1*/, int/*v0*/, int/*v1*/)>;

    /** \brief Check if the mesh is watertight (i.e., no holes or gaps). */
    using post_split_callback = std::function<void(const Eigen::MatrixXd& /*V*/, const Eigen::MatrixXi&/*F*/, 
                    const Eigen::MatrixXi &/*E*/, const Eigen::MatrixXi &/*uE*/, const Eigen::VectorXi& /*EMAP*/,
                    const std::vector<std::vector<int>> &/*uE2E*/, std::priority_queue<std::tuple<double,int,int>> &/*Q*/, 
                    int/*e0*/, int/*e1*/, int/*v0*/, int/*v1*/, int/*vM*/, int/*f0*/, int/*f1*/, int/*f2*/, int/*f3*/)>;

    /** \brief Check if the mesh is watertight (i.e., no holes or gaps). */
    bool splitter(const Eigen::MatrixXd& OV, const Eigen::MatrixXi& OF, split_condition_callback condition, 
                    split_placement_callback placement, pre_split_callback pre_split, post_split_callback post_split, 
                    Eigen::MatrixXd& V, Eigen::MatrixXi& F);

    /** \brief Check if the mesh is watertight (i.e., no holes or gaps). */
    bool is_watertight(const Eigen::MatrixXi& F);

    /** \brief Check if the mesh is watertight (i.e., no holes or gaps). */
    double compute_internal_angle(const Eigen::MatrixXd& V, int e0, int e1, int v);

    /** \brief Check if the mesh is watertight (i.e., no holes or gaps). */
    Eigen::RowVector3d project_vertex_to_edge(const Eigen::MatrixXd &V, int e0, int e1, int v_proj);

}