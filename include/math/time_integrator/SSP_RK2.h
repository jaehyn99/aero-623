#ifndef SSP_RK2_H
#define SSP_RK2_H

#include "TimeIntegrator.h"
class SSP_RK2: public TimeIntegrator{
    public:
    StateMesh integrate(const Function& f, const StateMesh& u0, double t, double dt) const override;
    StateMesh integrate(const Function& f, const StateMesh& u0, double t, const Eigen::ArrayXd& dt) const override;
};

#endif