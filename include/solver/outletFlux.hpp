#pragma once
#include "boundaryFlux.hpp"
#include <Eigen/Dense>

class outletFlux : public boundaryFlux {
private:
    double p0_;

public:
    outletFlux(double p0)
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
	
	Eigen::Matrix<double, 4, 1> F;

	if (std::sqrt(uB*uB + vB*vB)/cB < 1) {
		F(0) = rhoB*(uB*n(0) + vB*n(1));
		F(1) = (rhoB*uB*uB + p0_)*n(0) + rhoB*uB*vB*n(1);
		F(2) = rhoB*uB*vB*n(0) + (rhoB*vB*vB + p0_)*n(1);
		F(3) = uB*(rhoEB + p0_)*n(0) + vB*(rhoEB + p0_)*n(1);
	}
	else {
		F(0) = rhoP*(uP*n(0) + vP*n(1));
		F(1) = (rhoP*uP*uP + pP)*n(0) + rhoP*uP*vP*n(1);
		F(2) = rhoP*uP*vP*n(0) + (rhoP*vP*vP + pP)*n(1);
		F(3) = uP*(rhoEP + pP)*n(0) + vP*(rhoEP + pP)*n(1);
	}
	
        return F;
    }
};
