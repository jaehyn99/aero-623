#include "FreeStreamBC.h"

Eigen::Vector4d FreeStreamBC::computeFlux(const Eigen::Vector4d& U, const Eigen::Vector2d& n) const{
    double gm1 = _gamma - 1.0;
    double rho = U(0);
    double u = U(1)/rho;
    double v = U(2)/rho;
    double rhoE = U(3);
    double p = gm1*(rhoE - 0.5*rho*(u*u + v*v));
 
    Eigen::Matrix<double,4,2> F;

    F(0,0) = rho*u;
    F(1,0) = rho*u*u + p;
    F(2,0) = rho*u*v;
    F(3,0) = (rhoE + p)*u;

    F(0,1) = rho*v;
    F(1,1) = rho*u*v;
    F(2,1) = rho*v*v + p;
    F(3,1) = (rhoE + p)*v;
    return F.col(0)*n(0) + F.col(1)*n(1);
}