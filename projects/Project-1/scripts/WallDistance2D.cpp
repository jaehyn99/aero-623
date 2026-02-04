#include "mesh/WallDistance2D.h"

#include <limits>
#include <stdexcept>
#include <cmath>

namespace mesh {

std::vector<double> wallDistance2D(
    const std::vector<Eigen::Vector2d>& X,
    const std::vector<char>& isWall
){
    const std::size_t nNodes = X.size();
    if (isWall.size() != nNodes) {
        throw std::invalid_argument("wallDistance2D: isWall.size() != X.size()");
    }

    // Collect wall node coordinates
    std::vector<Eigen::Vector2d> wallPts;
    wallPts.reserve(nNodes);
    for (std::size_t i = 0; i < nNodes; ++i) {
        if (isWall[i]) wallPts.push_back(X[i]);
    }

    // For robustness: if no wall nodes exist, return +inf distances
    if (wallPts.empty()) {
        std::vector<double> df(nNodes, std::numeric_limits<double>::infinity());
        return df;
    }

    // Compute min distance to wallPts for each node
    std::vector<double> df(nNodes, std::numeric_limits<double>::infinity());
    for (std::size_t i = 0; i < nNodes; ++i) {
        double best2 = std::numeric_limits<double>::infinity();
        const Eigen::Vector2d& xi = X[i];
        for (const auto& xw : wallPts) {
            const double d2 = (xi - xw).squaredNorm();
            if (d2 < best2) best2 = d2;
        }
        df[i] = std::sqrt(best2);
    }
    return df;
}

} // namespace mesh
