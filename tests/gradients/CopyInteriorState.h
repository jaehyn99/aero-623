#pragma once

#include <Eigen/Dense>
#include "solver/boundaryState.hpp"

class CopyInteriorState : public boundaryState
{
public:
    Eigen::Vector4d operator()(
        const Eigen::Vector4d& UP,
        double,
        const Eigen::Vector2d&
    ) const override
    {
        return UP;   // critical
    }
};
