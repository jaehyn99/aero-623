#include "RK4.h"
#include "StateMesh.h"

void RK4::integrate(const Function& f, Eigen::MatrixXd& u0, double t, double dt) const{
    Eigen::MatrixXd k1 = dt*f(t, u0);

    Eigen::MatrixXd u_temp = u0 + 0.5*k1;
    Eigen::MatrixXd k2 = dt*f(t+0.5*dt, u_temp);

    u_temp = u0 + 0.5*k2;
    Eigen::MatrixXd k3 = dt*f(t+0.5*dt, u_temp);

    u_temp = u0 + k3;
    Eigen::MatrixXd k4 = dt*f(t+dt, u_temp);
    
    u0 = u0 + (k1+2*k2+2*k3+k4)/6;
}

void RK4::integrate(const Function& f, Eigen::MatrixXd& u0, double t, const Eigen::ArrayXd& dt) const{
    double dtmin = dt.minCoeff();
    Eigen::MatrixXd k1 = f(t, u0).array().rowwise() * dt.transpose();

    Eigen::MatrixXd u_temp = u0 + 0.5*k1;
    Eigen::MatrixXd k2 = f(t+0.5*dtmin, u_temp).array().rowwise() * dt.transpose();

    u_temp = u0 + 0.5*k2;
    Eigen::MatrixXd k3 = f(t+0.5*dtmin, u_temp).array().rowwise() * dt.transpose();

    u_temp = u0 + k3;
    Eigen::MatrixXd k4 = f(t+dtmin, u_temp).array().rowwise() * dt.transpose();

    u0 = u0 + (k1+2*k2+2*k3+k4)/6;
} 