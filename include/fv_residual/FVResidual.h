#ifndef FV_RESIDUAL_H
#define FV_RESIDUAL_H

#include "Eigen/Dense"
class StateMesh;
class FVResidual{
    public:
    virtual ~FVResidual() = default;
    virtual Eigen::MatrixXd computeResidual(const StateMesh& u) const = 0;
};

#endif