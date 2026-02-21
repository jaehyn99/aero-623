#pragma once

#include <Eigen/Dense>
#include "boundary_condition/BoundaryCondition.h"

class CopyBC : public BoundaryCondition {
public:
    Eigen::Vector4d computeFlux(const Eigen::Vector4d&, const Eigen::Vector2d&) const override {
        return Eigen::Vector4d::Zero(); // not used in gradient test
    }
    Eigen::Vector4d computeBoundaryState(const Eigen::Vector4d& UP, const Eigen::Vector2d&) const override {
        return UP; // consistent with any manufactured interior value
    }
};
