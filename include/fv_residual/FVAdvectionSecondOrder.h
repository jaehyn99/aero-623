#ifndef FV_ADVECTION_SECOND_ORDER_H
#define FV_ADVECTION_SECOND_ORDER_H

#include "FVResidual.h"

using GradMatrix = std::vector<Eigen::Matrix<double,4,2>>;

class FVFlux;
class FVAdvectionSecondOrder: public FVResidual{
    public:
    FVAdvectionSecondOrder(std::shared_ptr<FVFlux> flux): _flux(flux) {}
    Eigen::MatrixXd computeResidual(const StateMesh& u) const override;

    protected:
    std::shared_ptr<FVFlux> _flux;
};

#endif