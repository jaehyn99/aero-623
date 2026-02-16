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

        // Safe Guard
        constexpr double pFloor = 1e-10;
        pL = std::max(pFloor, pL);
        pR = std::max(pFloor, pR);


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

        double aL = std::sqrt(gamma * pL / rhoL);
        double aR = std::sqrt(gamma * pR / rhoR);

        double SL = std::min(uNL - aL, uNR - aR);
        double SR = std::max(uNL + aL, uNR + aR);

        if (SL >= 0.0) return FL;
        if (SR <= 0.0) return FR;
        
        F = (SR*FL - SL*FR + SL*SR*(UR - UL)) / (SR - SL);
        return F;
    }
};
