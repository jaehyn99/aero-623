#ifndef FV_ADVECTION_SECOND_ORDER_H
#define FV_ADVECTION_SECOND_ORDER_H

#include "FVResidual.h"

using GradMatrix = std::vector<Eigen::Matrix<double,4,2>>;

class FVFlux;
class FVAdvectionSecondOrder: public FVResidual{
    public:
    FVAdvectionSecondOrder(std::shared_ptr<FVFlux> flux, bool limiter=true): _flux(flux), _limiter(limiter) {}
    Eigen::MatrixXd computeResidual(const StateMesh& u) const override { return Eigen::MatrixXd::Zero(1,1); };

    protected:
    std::shared_ptr<FVFlux> _flux;
    bool _limiter;
};

#endif