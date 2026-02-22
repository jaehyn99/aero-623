#ifndef FV_GRADIENT_H
#define FV_GRADIENT_H

#include "Eigen/Dense"
class StateMesh;
class FVGradient
{
public:

    // Plane-normal gradient construction with hybrid treatment:
    // - Interior and boundary faces: walk/Green-Gauss style accumulation using In

    virtual std::vector<Eigen::Matrix<double,4,2>> computeGradient(const StateMesh& U) const = 0;
};

#endif