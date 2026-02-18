#include "SSP_RK2.h"
#include <iostream>

void SSP_RK2::integrate(const Function& f, Eigen::MatrixXd& u0, double t, double dt) const{
    Eigen::MatrixXd u1 = u0 + dt*f(t, u0);
    Eigen::MatrixXd u2 = 0.5*u0 + 0.5*u1 + 0.5*dt*f(t+dt, u1);
    u0 = std::move(u2);
}

void SSP_RK2::integrate(const Function& f, Eigen::MatrixXd& u0, double t, const Eigen::ArrayXd& dt) const{
    double dtmin = dt.minCoeff();
    Eigen::MatrixXd u1 = u0 + (f(t, u0).array().rowwise() * dt.transpose()).matrix();
    Eigen::MatrixXd u2 = 0.5*u0 + 0.5*u1 + 0.5*(f(t+dtmin, u1).array().rowwise() * dt.transpose()).matrix();
    u0 = std::move(u2);
}