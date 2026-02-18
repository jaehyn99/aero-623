#include "SSP_RK2.h"
#include "StateMesh.h"
#include <iostream>

void SSP_RK2::integrate(const Function& f, StateMesh& u0, double t, double dt) const{
    Eigen::MatrixXd u1 = u0.matrix() + dt*f(t, u0.matrix());
    Eigen::MatrixXd u2 = 0.5*u0.matrix() + 0.5*u1 + 0.5*dt*f(t+dt, u1);

    // for (Eigen::Index i = 0; i < u2.cols(); i+=4){
    //     if (u2.col(i).array().isNaN().any()){
    //         auto mesh = u0.mesh();
    //         std::cout << u0.cell(i).transpose() << std::endl;
    //         std::cout << u1.col(i).transpose() << std::endl;
    //         break;
    //     }
    // }

    u0.matrix() = std::move(u2);
}

void SSP_RK2::integrate(const Function& f, StateMesh& u0, double t, const Eigen::ArrayXd& dt) const{
    double dtmin = dt.minCoeff();
    Eigen::MatrixXd u1 = u0.matrix() + (f(t, u0.matrix()).array().colwise() * dt).matrix();
    Eigen::MatrixXd u2 = 0.5*u0.matrix() + 0.5*u1 + 0.5*(f(t+dtmin, u1).array().colwise() * dt).matrix();
    u0.matrix() = std::move(u2);
}