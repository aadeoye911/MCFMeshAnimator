
#include "utils.h"

using namespace Eigen;

using ArrayXb = Eigen::Array<bool, Eigen::Dynamic,1>;

namespace utils {

bool splitter(const MatrixXd& OV, const MatrixXi& OF, split_condition_callback condition, split_placement_callback placement, 
                pre_split_callback pre_split, post_split_callback post_split, MatrixXd& V, MatrixXi& F) {

    V = OV;
    F = OF;

    MatrixXi E, uE;
    VectorXi EMAP;
    std::vector<std::vector<int>> uE2E;
    igl::unique_edge_map(F, E, uE, EMAP, uE2E);
    std::priority_queue<std::tuple<double,int,int>> Q;

    int vertex_count = V.rows();
    int face_count = F.rows();
    int edge_count = uE.rows();

    // Variables to track preallocated spaces
    int vertex_capacity = vertex_count;
    int face_capacity = face_count;
    int edge_capacity = edge_count;

    // Lambda function to preallocate memory in batches
    auto preallocate_memory = [&]() -> void {
        int queue_size = Q.size();
        // Update uE2E entries to align with updated indexing
        for (int j = 0; j < edge_count; j++) {
            uE2E[j][0] +=  (uE2E[j][0] / face_count) * 2 * queue_size;
            uE2E[j][1] +=  (uE2E[j][1] / face_count) * 2 * queue_size;
        }

        vertex_capacity = vertex_count + queue_size;
        face_capacity = face_count + 2 * queue_size;
        edge_capacity = edge_count + 3 * queue_size;
  
        V.conservativeResize(vertex_capacity, 3); 
        V.bottomRows(queue_size).setConstant(-1); 

        F.conservativeResize(face_capacity, 3);
        F.bottomRows(2*queue_size).setConstant(-1); 

        uE.conservativeResize(edge_capacity, 2);
        uE2E.resize(edge_capacity, std::vector<int>(2));

        // Update E and EMAP for new halfedge indexing
        int face_capacity = F.rows();
        MatrixXi E_temp(2 * edge_capacity, 2);
        VectorXi EMAP_temp(2 * edge_capacity);
        for (int i = 0; i < 3; i++) {
            E_temp.block(i * face_capacity, 0, face_count, 2) = E.block(i * face_count, 0, face_count, 2);
            E_temp.block(i * face_capacity + face_count, 0, 2 * queue_size, 2).setConstant(-1);  // Move block
            EMAP_temp.segment(i * face_capacity, face_count) = EMAP.segment(i * face_count, face_count);    // Move corresponding EMAP block
            EMAP_temp.segment(i * face_capacity + face_count,  2 * queue_size).setConstant(- 1);  // Set to last valid edge
        }
        E = E_temp;
        EMAP = EMAP_temp;
    };

    // Load edges to queue with conditional filter
    std::tuple<double,int,int> p;
    for (int uei = 0; uei < uE.rows(); uei++) {
        condition(V, F, E, uE, EMAP, uE2E, Q, uei);
    }
    std::cout << "Initial queue size: " << Q.size() << std::endl;

    while(!Q.empty()) {
        if (vertex_count == vertex_capacity) {
            preallocate_memory();
        }

        p = Q.top();
        Q.pop();
        int uei = std::get<1>(p);
        
        assert (uE2E[uei].size() == 2);
        int f0 = uE2E[uei][0] % face_capacity;
        int f1 = uE2E[uei][1] % face_capacity;
        int c0 = uE2E[uei][0] / face_capacity;
        int c1 = uE2E[uei][1] / face_capacity;
        assert(f0 != f1 && f0 != -1 && f1 != -1);
        assert(c0 < 3 && c1 < 3);

        int e0 = F(f0, (c0+1)%3); // first vertex of split edge
        int e1 = F(f0, (c0+2)%3); // second vertex of split edge
        int v0 = F(f0, c0);
        int v1 = F(f1, c1);
        assert(F(f1, (c1+1)%3) == e1);
        assert(F(f1, (c1+2)%3) == e0);
        assert(EMAP(uE2E[uei][0]) == uei);
        assert(EMAP(uE2E[uei][1]) == uei);

        // Callback to update the queue if dynamic condition
        if (!pre_split(V, F, E, uE, EMAP, uE2E, Q, p, e0, e1, v0, v1)) {
            continue;
        }

        int vM = vertex_count; // midpoint vertex
        V.row(vM) = placement(V, F, E, uE, EMAP, uE2E, p, e0, e1, v0, v1);

        // *** UPDATE F ***
        int f2 = face_count;
        int f3 = face_count + 1;
        F.row(f2) = F.row(f0);                         // [e0, e1, v0] in some cyclic order
        F.row(f3) = F.row(f1);                         // [e1, e0, v1] in some cyclic order
        // Replacement subfaces adjacent to [e0, vM]
        F(f0, (c0+2) % 3) = vM;                        // f0: [e0, e1, v0] -> [e0, vM, v0]
        F(f1, (c1+1) % 3) = vM;                        // f1: [e1, e0, v1] -> [vM, e0, v1]
        // Additional subfaces adjacent to [vM, e1]
        F(f2, (c0+1) % 3) = vM;                        // f2: [e0, e1, v0] -> [vM, e1, v0]
        F(f3, (c1+2) % 3) = vM;                        // f3: [e1, e0, v1] -> [e1, vM, v1]

        // *** UPDATE uE ***
        uE.row(uei) << e0, vM; // Replacement edge [e0, e1] -> [e0, vM]
        uE.row(edge_count) << vM, v0;
        uE.row(edge_count + 1) << vM, v1;
        uE.row(edge_count + 2) << vM, e1;

        // *** UPDATE E ***
        E.row(f0 + face_capacity * c0) << e0, vM;
        E.row(f1 + face_capacity * c1) << vM, e0;
        E.row(f0 + face_capacity * ((c0+1) % 3)) << vM, v0;
        E.row(f1 + face_capacity * ((c1+2) % 3)) << v1, vM;
        E.row(f2 + face_capacity * ((c0+1) % 3)) << e1, v0;
        E.row(f3 + face_capacity * ((c1+2) % 3)) << v1, e1;
        E.row(f2 + face_capacity * ((c0+2) % 3)) << v0, vM;
        E.row(f3 + face_capacity * ((c1+1) % 3)) << vM, v1;
        E.row(f2 + face_capacity * c0) << vM, e1;
        E.row(f3 + face_capacity * c1) << e1, vM;

        // *** UPDATE EMAP ***
        int uei_e1v0 = EMAP(f0 + face_capacity * ((c0+1) % 3));
        int uei_v1e1 = EMAP(f1 + face_capacity * ((c1+2) % 3));
        EMAP(f0 + face_capacity * ((c0+1) % 3)) = edge_count;     // [vM, v0]
        EMAP(f1 + face_capacity * ((c1+2) % 3)) = edge_count + 1; // [v1, vM]
        EMAP(f2 + face_capacity * c0) = edge_count + 2;           // [vM, e1]
        EMAP(f2 + face_capacity * ((c0+1) % 3)) = uei_e1v0;       // [e1, v0]
        EMAP(f2 + face_capacity * ((c0+2) % 3)) = edge_count;     // [v0, vM]
        EMAP(f3 + face_capacity * c1) = edge_count + 2;           // [e1, vM]
        EMAP(f3 + face_capacity * ((c1+1) % 3)) = edge_count + 1; // [vM, v1]
        EMAP(f3 + face_capacity * ((c1+2) % 3)) = uei_v1e1;       // [v1, e1]

        int uei_v0e0 = EMAP(f0 + face_capacity * ((c0+2) % 3));
        int uei_e0v1 = EMAP(f1 + face_capacity * ((c1+1) % 3));
        assert(uE2E[uei_v0e0][0] == f0 + face_capacity * ((c0+2) % 3) || uE2E[uei_v0e0][1] == f0 + face_capacity * ((c0+2) % 3));
        assert(uE2E[uei_e0v1][0] == f1 + face_capacity * ((c1+1) % 3) || uE2E[uei_e0v1][1] == f1 + face_capacity * ((c1+1) % 3));
        
        // *** UPDATE uE2E ***
        uE2E[edge_count][0] = f0 + face_capacity * ((c0+1) % 3);
        uE2E[edge_count][1] = f2 + face_capacity * ((c0+2) % 3);
        uE2E[edge_count + 1][0] = f1 + face_capacity * ((c1+2) % 3);
        uE2E[edge_count + 1][1] = f3 + face_capacity * ((c1+1) % 3);
        uE2E[edge_count + 2][0] = f2 + face_capacity * c0;
        uE2E[edge_count + 2][1] = f3 + face_capacity * c1;

        if (uE2E[uei_e1v0][0] == f0 + face_capacity * ((c0+1) % 3)) {
            uE2E[uei_e1v0][0] = f2 + face_capacity * ((c0+1) % 3);
        }
        else {
            assert(uE2E[uei_e1v0][1] == f0 + face_capacity * ((c0+1) % 3));
            uE2E[uei_e1v0][1] = f2 + face_capacity * ((c0+1) % 3);
        }

        if (uE2E[uei_v1e1][0] == f1 + face_capacity * ((c1+2) % 3)) {
            uE2E[uei_v1e1][0] = f3 + face_capacity * ((c1+2) % 3);
        }
        else {
            assert(uE2E[uei_v1e1][1] == f1 + face_capacity * ((c1+2) % 3));
            uE2E[uei_v1e1][1] = f3 + face_capacity * ((c1+2) % 3);
        }

        // Callback to update the queue follow edge splits
        post_split(V, F, E, uE, EMAP, uE2E, Q, e0, e1, v0, v1, vM, f0, f1, f2, f3);

        // Update counters
        vertex_count++;
        face_count += 2;
        edge_count += 3;
    }

    if (V.rows() != vertex_count) {
        assert(V.rows() > vertex_count);
        assert(vertex_count == F.maxCoeff() + 1);
        assert((V.bottomRows(V.rows() - vertex_count - 1).array() == -1).all());
        assert((F.bottomRows(F.rows() - face_count).array() == -1).all());
        
        // Remove the excess rows from V and F
        V.conservativeResize(vertex_count, V.cols());  
        F.conservativeResize(face_count, F.cols()); 
    }

    return utils::is_watertight(F);
}

double compute_internal_angle(const MatrixXd& V, int e0, int e1, int w) {
    auto A = (V.row(e0) - V.row(w)).normalized();
    auto B = (V.row(e1) - V.row(w)).normalized(); 
    
    return 2.0 * atan((A - B).norm() / (A + B).norm());
};

RowVector3d project_vertex_to_edge(const MatrixXd& V, int e0, int e1, int v) {
    RowVector3d p0 = V.row(e0);  // Point e0
    RowVector3d p1 = V.row(e1);  // Point e1
    RowVector3d p = V.row(v);   // Point v0 to project onto the edge

    RowVector3d edge = p1 - p0;
    RowVector3d projectee = p - p0;

    double t = projectee.dot(edge) / edge.squaredNorm();  // Compute the projection factor
    t = std::max(0.0, std::min(1.0, t));
    
    return p0 + t * edge;
}

bool is_watertight(const MatrixXi& F) {
    if (!igl::is_vertex_manifold(F)) {
        std::cerr << "Mesh is not vertex manifold" << std::endl;
        return false;
    }
    if (!igl::is_edge_manifold(F)) {
        std::cerr << "Mesh is not edge manifold" << std::endl;
        return false;
    }
    std::vector<std::vector<int>> boundary_loops;
    igl::boundary_loop(F, boundary_loops);
    if (!boundary_loops.empty()) {
        std::cerr << "Mesh has an open boundary" << std::endl;
        return false;
    }
    return true;
}

} // Namespace 