#pragma once

#include <Eigen/Dense>
#include "solver/boundaryState.hpp"

// --------------------------------------------------
// Constant boundary state for testing
// Returns U = [1,1,1,1]^T regardless of input
// --------------------------------------------------

class ConstantBoundaryState : public boundaryState
{
public:
    ConstantBoundaryState() = default;

    Eigen::Vector4d operator()(
        const Eigen::Vector4d& /*UP*/,
        double /*gamma*/,
        const Eigen::Vector2d& /*n*/
    ) const override
    {
        return Eigen::Vector4d::Ones();
    }
};
