#pragma once

#include <Eigen/Dense>
#include <vector>

#include "mesh/GriIO.h"
#include "mesh/Projection2D.h"   // Projection2D + ProjectionResult
#include "mesh/BladeGeometry.h"  // BladeGeometry

namespace mesh {

class SizingFunction; // forward declare is OK here

class LocalRefinement {
public:
    static void computeNodeDistanceAndSize(
        const Eigen::MatrixXd& V,
        const BladeGeometry& blade,
        const Projection2D& projector,
        const SizingFunction& sizeFun,
        Eigen::VectorXd& distOut,
        Eigen::VectorXd& hOut
    );

    static ProjectionResult projectDistance(
        const Eigen::Vector2d& p,
        const Eigen::MatrixXd& V,   // only used to get Ly
        const BladeGeometry& blade,
        const Projection2D& projector
    );

    static void markElementsForRefinement(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::VectorXd& hNode,
        double alpha,
        std::vector<char>& refineFlag
    );

    static GriMesh2D refineMesh(
        const GriMesh2D& meshIn,
        const BladeGeometry& blade,
        const Projection2D& projector,
        const SizingFunction& sizeFun,
        double alpha,
        int maxIters
    );

    static int refineAdaptive(
        Eigen::MatrixXd& V,
        Eigen::MatrixXi& E,
        const Eigen::MatrixXi& boundaryEdges,
        const BladeGeometry& blade,
        const Projection2D& projector,
        const SizingFunction& sizeFun,
        double C = 1.25,
        int maxIters = 20,
        int smoothIters = 5,
        double omega = 0.5
    );
};

} // namespace mesh
