#include "SSP_RK3.h"
#include "StateMesh.h"

void SSP_RK3::integrate(const Function& f, Eigen::MatrixXd& u0, double t, double dt) const{
    Eigen::MatrixXd u1 = u0 + dt*f(t, u0);
    Eigen::MatrixXd u2 = 0.75*u0 + 0.25*u1 + 0.25*dt*f(t+dt, u1);
    Eigen::MatrixXd u3 = u0/3 + 2*u2/3 + 2*dt*f(t+0.5*dt, u2)/3;
    u0 = std::move(u3);
}

void SSP_RK3::integrate(const Function& f, Eigen::MatrixXd& u0, double t, const Eigen::ArrayXd& dt) const{
    double dtmin = dt.minCoeff();
    Eigen::MatrixXd u1 = u0 + (f(t, u0).array().rowwise() * dt.transpose()).matrix();
    Eigen::MatrixXd u2 = 3*u0/4 + u1/4 + ((f(t+dtmin, u1)/4).array().rowwise() * dt.transpose()).matrix();
    Eigen::MatrixXd u3 = u0/3 + 2*u2/3 + 2.0/3*(f(t+0.5*dtmin, u2).array().rowwise() * dt.transpose()).matrix();
    u0 = std::move(u3);
}