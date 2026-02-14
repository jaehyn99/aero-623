#pragma once
#include "boundaryState.hpp"
#include <Eigen/Dense>

class outletState : public boundaryState {
private:
    double p0_;

public:
    outletState(double p0)
        : p0_(p0) {}

    Eigen::Vector4d operator()(
        const Eigen::Vector4d& UP,
        double gamma,
        const Eigen::Vector2d& n
    ) const override {
	double gm1 = gamma - 1;
        double rhoP = UP(0); 
	double uP = UP(1)/rhoP;
	double vP = UP(2)/rhoP;
	double rhoEP = UP(3);
	double pP = gm1*(rhoEP - 0.5*rhoP*(uP*uP + vP*vP));
	double SP = p0_/std::pow(rhoP, gamma);
	double rhoB = std::pow(p0_/SP, 1/gamma);
	double cB = std::sqrt(gamma*p0_/rhoB);
	double uNP = uP*n(0) + vP*n(1);
	double cP = std::sqrt(gamma*pP/rhoP);
	double JP = uNP + 2*cP/gm1;
	double uNB = JP - 2*cB/gm1;
	double uB = uP - n(0)*(uP*n(0) + vP*n(1) - uNB);
	double vB = vP - n(1)*(uP*n(0) + vP*n(1) - uNB);
	double rhoEB = p0_/gm1 + 0.5*rhoB*(uB*uB + vB*vB);
	
	Eigen::Matrix<double, 4, 1> Ub;

	if (std::sqrt(uB*uB + vB*vB)/cB < 1) {
		Ub(0) = rhoB;
		Ub(1) = rhoB*uB;
		Ub(2) = rhoB*vB;
		Ub(3) = rhoEB;
	}
	else {
		Ub(0) = rhoP;
		Ub(1) = rhoP*uP;
		Ub(2) = rhoP*vP;
		Ub(3) = rhoEP;
	}
	
        return Ub;
    }
};
