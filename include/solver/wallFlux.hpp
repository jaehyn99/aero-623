#pragma once
#include "boundaryFlux.hpp"
#include <Eigen/Dense>

class wallFlux : public boundaryFlux {
public:
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
	double uB = uP - n(0)*(uP*n(0) + vP*n(1));
	double vB = vP - n(1)*(uP*n(0) + vP*n(1));
	double pB = gm1*(rhoEP - 0.5*rhoP*(uB*uB + vB*vB));
	double cB = std::sqrt(gamma*pB/rhoP);
	double pP = gm1*(rhoEP - 0.5*rhoP*(uP*uP + vP*vP));
	double cP = std::sqrt(gamma*pP/rhoP);
	double uNP = uP*n(0) + vP*n(1);
	double JP = uNP + 2*cP/gm1;
	double uNB = JP - 2*cB/gm1;
	Eigen::Matrix<double, 4, 1> F;
	F(0) = 0.0;
	F(1) = pB*n(0);
	F(2) = pB*n(1);
	F(3) = 0.0;
        return F;
    }
};
