#pragma once
#include <Eigen/Dense>
#include <string>
#include <vector>

namespace mesh {

struct GriBoundaryBlock {
    int tag = 0;
    std::string name;
    Eigen::MatrixXi edges; // Nb x 2 (0-based)
};

struct GriMesh2D {
    Eigen::MatrixXd V;                 // Nn x 2
    Eigen::MatrixXi E;                 // Ne x 3
    std::vector<GriBoundaryBlock> B;   // boundary blocks
};

GriMesh2D readGri2D(const std::string& fname);
void writeGri2D(const std::string& fname, const GriMesh2D& mesh);

} // namespace mesh
