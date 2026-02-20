#include "HLLEFlux.h"
#include <iostream>

Eigen::Vector4d HLLEFlux::computeFlux(const Eigen::Vector4d& UL, const Eigen::Vector4d& UR, const Eigen::Vector2d& n) const{
	double gm1 = _gamma - 1.0;

    double rhoL = UL(0);
    double uL = UL(1)/rhoL;
    double vL = UL(2)/rhoL;
    double uNL = uL*n(0) + vL*n(1);
    double pL = gm1*(UL(3) - 0.5*rhoL*(uL*uL + vL*vL));

    double rhoR = UR(0);
    double uR = UR(1)/rhoR;
    double vR = UR(2)/rhoR;
	double uNR = uR*n(0) + vR*n(1);    
    double pR = gm1*(UR(3) - 0.5*rhoR*(uR*uR + vR*vR));

	Eigen::Vector4d FL{rhoL*uNL, UL(1)*uNL + pL*n(0), UL(2)*uNL + pL*n(1), (UL(3) + pL)*uNL};
    Eigen::Vector4d FR{rhoR*uNR, UR(1)*uNR + pR*n(0), UR(2)*uNR + pR*n(1), (UR(3) + pR)*uNR};

    double eta2 = 0.5*std::sqrt(rhoL*rhoR) / std::pow(std::sqrt(rhoL) + std::sqrt(rhoR), 2);
    double vAvg = (std::sqrt(rhoL)*uNL + std::sqrt(rhoR)*uNR) / (std::sqrt(rhoL) + std::sqrt(rhoR));
    double EL = UL(3)/rhoL;
    double ER = UR(3)/rhoR;

    double aL2 = gm1*((EL + pL)/rhoL - 0.5*(uL*uL + vL*vL));
    double aR2 = gm1*((ER + pR)/rhoR - 0.5*(uR*uR + vR*vR));

    double dAvg = (std::sqrt(rhoL)*aL2 + std::sqrt(rhoR)*aR2) / (std::sqrt(rhoL) + std::sqrt(rhoR)) + eta2*(uNR - uNL)*(uNR - uNL);
    //std::cout << dAvg << ", " << vAvg << std::endl;
    dAvg = std::sqrt(dAvg);
    double SL = vAvg - dAvg;
    double SR = vAvg + dAvg;

    if (SL >= 0.0) return FL; // Want to do masking later
    if (SR <= 0.0) return FR;
    return (SR*FL - SL*FR + SL*SR*(UR-UL))/(SR-SL);
}

double HLLEFlux::computeWaveSpeed(const Eigen::Vector4d& UL, const Eigen::Vector4d& UR, const Eigen::Vector2d& n, double cL, double cR) const{
    double uL = std::abs(UL(1)/UL(0)*n(0) + UL(2)/UL(0)*n(1));
    double uR = std::abs(UR(1)/UR(0)*n(0) + UR(2)/UR(0)*n(1));
    return std::max(uL+cL, uR+cR);
}