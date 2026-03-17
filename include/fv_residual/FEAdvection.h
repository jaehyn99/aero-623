#ifndef FE_ADVECTION_H
#define FE_ADVECTION_H

#include "FVResidual.h"

class FVFlux;
class FEAdvection: public FVResidual{
    public:
    FEAdvection(std::shared_ptr<FVFlux> flux): _flux(flux) {}
    Eigen::MatrixXd computeResidual(const StateMesh& u) const override;

    protected:
    std::shared_ptr<FVFlux> _flux;
};

#endif