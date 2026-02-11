#pragma once
#include <Eigen/Dense>

class Flux {
public:
    virtual Eigen::Vector4d operator()(const Eigen::Vector4d& UL,
                     const Eigen::Vector4d& UR,
                     double gamma, const Eigen::Vector2d& n) const = 0;
    
    virtual ~Flux() = default;
};
