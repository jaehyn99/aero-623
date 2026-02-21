#pragma once
#include <Eigen/Dense>
#include "boundary_condition/BoundaryCondition.h"

// Constant boundary condition for the BoundaryCondition/StateMesh framework
class ConstantBC : public BoundaryCondition {
public:
    explicit ConstantBC(const Eigen::Vector4d& Uconst = Eigen::Vector4d::Ones())
        : Uc_(Uconst) {}

    Eigen::Vector4d computeFlux(const Eigen::Vector4d&, const Eigen::Vector2d&) const override {
        return Eigen::Vector4d::Zero(); // not needed for gradient tests
    }

    Eigen::Vector4d computeBoundaryState(const Eigen::Vector4d&, const Eigen::Vector2d&) const override {
        return Uc_;
    }

private:
    Eigen::Vector4d Uc_;
};
