#ifndef SSP_RK2_H
#define SSP_RK2_H

#include "TimeIntegrator.h"
class SSP_RK2: public TimeIntegrator{
    public:
    void integrate(const Function& f, Eigen::MatrixXd& u0, double t, double dt) const override;
    void integrate(const Function& f, Eigen::MatrixXd& u0, double t, const Eigen::ArrayXd& dt) const override;
};

#endif