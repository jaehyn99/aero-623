#ifndef SSP_RK3_H
#define SSP_RK3_H

#include "TimeIntegrator.h"
class SSP_RK3: public TimeIntegrator{
    public:
    void integrate(const Function& f, StateMesh& u0, double t, double dt) const override;
    void integrate(const Function& f, StateMesh& u0, double t, const Eigen::ArrayXd& dt) const override;
};

#endif