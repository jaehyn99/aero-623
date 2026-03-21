#ifndef RK4_H
#define RK4_H

#include "TimeIntegrator.h"
class RK4: public TimeIntegrator{
    public:
    void integrate(const Function& f, Eigen::MatrixXd& u0, double t, double dt) const override;
    void integrate(const Function& f, Eigen::MatrixXd& u0, double t, const Eigen::ArrayXd& dt) const override;
};

#endif