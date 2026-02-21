#pragma once
#include <Eigen/Dense>

class boundaryFlux {
public:
    virtual FluxResult operator()(const Eigen::Vector4d& UP,
                     double gamma, const Eigen::Vector2d& n) const = 0;
    
    virtual ~boundaryFlux() = default;
};
