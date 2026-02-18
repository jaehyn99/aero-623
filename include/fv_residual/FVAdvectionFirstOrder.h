#ifndef FV_ADVECTION_FIRST_ORDER_H
#define FV_ADVECTION_FIRST_ORDER_H

#include "FVResidual.h"

class FVFlux;
class FVAdvectionFirstOrder: public FVResidual{
    public:
    FVAdvectionFirstOrder(std::shared_ptr<FVFlux> flux): _flux(flux) {}
    Eigen::MatrixXd computeResidual(const StateMesh& u) const override;

    protected:
    std::shared_ptr<FVFlux> _flux;
};

#endif