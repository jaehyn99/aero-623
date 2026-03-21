#ifndef HYBRID_WALK_PN_GRAD_H
#define HYBRID_WALK_PN_GRAD_H

#include "FVGradient.h"

class HybridWalkPNGrad : public FVGradient
{
public:

    // Plane-normal gradient construction with hybrid treatment:
    // - Interior faces: walk/Green-Gauss style accumulation using In
    // - Boundary faces: plane-normal reconstruction using two interior neighbors
    std::vector<Eigen::Matrix<double,4,2>> computeGradient(const StateMesh& U) const override{ return {}; };
};

#endif