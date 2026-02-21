#pragma once
#include <Eigen/Dense>
#include "FluxResult.hpp"

class numericalFlux {
public:
    virtual FluxResult operator()(const Eigen::Vector4d& UL,
                     const Eigen::Vector4d& UR,
                     double gamma, const Eigen::Vector2d& n) const = 0;
    
    virtual ~numericalFlux() = default;
};
