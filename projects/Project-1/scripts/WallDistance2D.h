#pragma once

#include <vector>
#include <Eigen/Dense>

namespace mesh {

/**
 * Compute nodal wall distance in 2D as the minimum Euclidean distance
 * to the set of wall nodes (nearest wall-node distance).
 *
 * Inputs:
 *  - X: nodal coordinates (size nNodes)
 *  - isWall: marker array (size nNodes), 1 if node is a wall node
 *
 * Output:
 *  - df: wall distance per node (size nNodes)
 *
 * Notes:
 *  - This matches the legacy Fortran "walldist" behavior (node-to-wall-node).
 *  - If you want exact wall distance to boundary edges, use a segment-based method.
 */
std::vector<double> wallDistance2D(
    const std::vector<Eigen::Vector2d>& X,
    const std::vector<char>& isWall
);

} // namespace mesh
