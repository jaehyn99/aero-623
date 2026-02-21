#pragma once

#include "mesh/TriangularMesh.h"
#include "mesh/StateMesh.h"
#include "Gradients.h"
#include <Eigen/Dense>

class HybridWalkPNGrad : public Gradients
{
public:

    // Plane-normal gradient construction with hybrid treatment:
    // - Interior faces: walk/Green-Gauss style accumulation using In
    // - Boundary faces: plane-normal reconstruction using two interior neighbors
    std::vector<Eigen::Matrix<double,4,2>> computeGradient(const Eigen::MatrixXi& I2E,
                                                        const Eigen::MatrixXi& B2E,
                                                        const Eigen::MatrixXd& In,
                                                        const Eigen::MatrixXd& Bn,
                                                        const Eigen::VectorXd& Area,
                                                        const StateMesh& U) const override;
};
