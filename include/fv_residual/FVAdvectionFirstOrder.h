#ifndef FV_ADVECTION_FIRST_ORDER_H
#define FV_ADVECTION_FIRST_ORDER_H

#include "Residual.h"

class FVFlux;
class FVAdvectionFirstOrder: public Residual{
    public:
    FVAdvectionFirstOrder(std::shared_ptr<FVFlux> flux): _flux(flux) {}
    Eigen::MatrixXd computeResidual(const StateMesh& u) const override { return Eigen::MatrixXd::Zero(1,1); };

    protected:
    std::shared_ptr<FVFlux> _flux;
};

#endif