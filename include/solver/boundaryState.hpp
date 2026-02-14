#pragma once
#include <Eigen/Dense>

class boundaryState {
public:
    virtual Eigen::Vector4d operator()(const Eigen::Vector4d& UP,
                     double gamma, const Eigen::Vector2d& n) const = 0;
    
    virtual ~boundaryState() = default;
};
