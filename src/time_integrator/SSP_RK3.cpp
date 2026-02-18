#include "SSP_RK3.h"
#include "StateMesh.h"

void SSP_RK3::integrate(const Function& f, StateMesh& u0, double t, double dt) const{
    Eigen::MatrixXd u1 = u0.matrix() + dt*f(t, u0.matrix());
    Eigen::MatrixXd u2 = 0.75*u0.matrix() + 0.25*u1 + 0.25*dt*f(t+dt, u1);
    Eigen::MatrixXd u3 = u0.matrix()/3 + 2*u2/3 + 2*dt*f(t+0.5*dt, u2)/3;
    u0.matrix() = std::move(u3);
}

void SSP_RK3::integrate(const Function& f, StateMesh& u0, double t, const Eigen::ArrayXd& dt) const{
    double dtmin = dt.minCoeff();
    Eigen::MatrixXd u1 = u0.matrix() + (f(t, u0.matrix()).array().colwise() * dt).matrix();
    Eigen::MatrixXd u2 = 0.75*u0.matrix() + 0.25*u1 + 0.25*(f(t+dtmin, u1).array().colwise() * dt).matrix();
    Eigen::MatrixXd u3 = u0.matrix()/3 + 2*u2/3 + 2/3*(f(t+0.5*dtmin, u2).array().colwise() * dt).matrix();
    u0.matrix() = std::move(u3);
}