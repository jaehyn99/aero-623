#pragma once
#include "numericalFlux.hpp"
#include <Eigen/Dense>
#include <cmath>

// First-order solver-specific HLLE implementation (kept numerically identical
// to legacy HLLEFlux for controlled A/B checks).
class HLLEFluxFO : public numericalFlux {
public:
    Eigen::Vector4d operator()(const Eigen::Vector4d& UL,
                     const Eigen::Vector4d& UR,
                     double gamma, const Eigen::Vector2d& n) const override
    {
	double gm1 = gamma - 1.0;

        double rhoL = UL(0);
        double rhoR = UR(0);

        double uL = UL(1)/rhoL;
        double uR = UR(1)/rhoR;

        double vL = UL(2)/rhoL;
        double vR = UR(2)/rhoR;

	double uNL = uL*n(0) + vL*n(1);
	double uNR = uR*n(0) + vR*n(1);

        double pL = gm1*(UL(3) - 0.5*rhoL*(uL*uL + vL*vL));
        double pR = gm1*(UR(3) - 0.5*rhoR*(uR*uR + vR*vR));

	Eigen::Matrix<double, 4, 1> FL; Eigen::Matrix<double, 4, 1> FR; 
	Eigen::Matrix<double, 4, 1> F;

	FL(0) = rhoL*uNL;
	FL(1) = UL(1)*uNL + pL*n(0);
	FL(2) = UL(2)*uNL + pL*n(1);
	FL(3) = (UL(3) + pL)*uNL;

	FR(0) = rhoR*uNR;
	FR(1) = UR(1)*uNR + pR*n(0);
	FR(2) = UR(2)*uNR + pR*n(1);
	FR(3) = (UR(3) + pR)*uNR;

        double eta2 = 0.5*std::sqrt(rhoL*rhoR)
                      / std::pow(std::sqrt(rhoL) + std::sqrt(rhoR), 2);

        double vAvg = (std::sqrt(rhoL)*vL + std::sqrt(rhoR)*vR)
                      / (std::sqrt(rhoL) + std::sqrt(rhoR));

        double EL = UL(3)/rhoL;
        double ER = UR(3)/rhoR;

        double aL2 = gm1*((EL + pL)/rhoL - 0.5*(uL*uL + vL*vL));
        double aR2 = gm1*((ER + pR)/rhoR - 0.5*(uR*uR + vR*vR));

        double dAvg = std::sqrt(
            (std::sqrt(rhoL)*aL2 + std::sqrt(rhoR)*aR2)
            / (std::sqrt(rhoL) + std::sqrt(rhoR))
            + eta2*(vR - vL)*(vR - vL)
        );

        double SL = vAvg - dAvg;
        double SR = vAvg + dAvg;

        if (SL >= 0.0) return FL;
        if (SR <= 0.0) return FR;

	F(0) = (SR*FL(0) - SL*FR(0) + SL*SR*(UR(0) - UL(0)))/(SR - SL);
	F(1) = (SR*FL(1) - SL*FR(1) + SL*SR*(UR(1) - UL(1)))/(SR - SL);
	F(2) = (SR*FL(2) - SL*FR(2) + SL*SR*(UR(2) - UL(2)))/(SR - SL);
	F(3) = (SR*FL(3) - SL*FR(3) + SL*SR*(UR(3) - UL(3)))/(SR - SL);

        return F;
    }
};
