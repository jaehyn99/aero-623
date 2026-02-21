#ifndef BOUNDARY_CONDITION_H
#define BOUNDARY_CONDITION_H

#include "Eigen/Dense"
class BoundaryCondition{
    public:
    virtual ~BoundaryCondition() = default;
    // returns the Euler flux
    virtual Eigen::Vector4d computeFlux(const Eigen::Vector4d& UP, const Eigen::Vector2d& n) const = 0;
    // returns a vector of density, x-momentum, v-momentum, and speed of sound at the boundary
    virtual Eigen::Vector4d computeBoundaryState(const Eigen::Vector4d& UP, const Eigen::Vector2d& n) const = 0;
};

#endif