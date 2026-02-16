#pragma once
#include "boundaryFlux.hpp"
#include <Eigen/Dense>

class wallFlux : public boundaryFlux {
public:
	// Eigen::Vector4d operator()(const Eigen::Vector4d& UP,
	// 						double gamma,
	// 						const Eigen::Vector2d& n) const override {
	// 	constexpr double rhoFloor = 1e-10;
	// 	constexpr double pFloor   = 1e-10;

	// 	const double gm1  = gamma - 1.0;
	// 	const double rhoP = std::max(rhoFloor, UP(0));
	// 	const double uP   = UP(1) / rhoP;
	// 	const double vP   = UP(2) / rhoP;
	// 	const double rhoE = UP(3);

	// 	// remove normal velocity component (slip wall)
	// 	const double un = uP*n(0) + vP*n(1);
	// 	const double uB = uP - un*n(0);
	// 	const double vB = vP - un*n(1);

	// 	double pB = gm1 * (rhoE - 0.5*rhoP*(uB*uB + vB*vB));
	// 	pB = std::max(pFloor, pB);

	// 	Eigen::Vector4d F;
	// 	F(0) = 0.0;
	// 	F(1) = pB * n(0);
	// 	F(2) = pB * n(1);
	// 	F(3) = 0.0;
	// 	return F;
	// }

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
				
	constexpr double pFloor = 1e-10;
	double pB = gm1*(rhoEP - 0.5*rhoP*(uB*uB + vB*vB));
	pB = std::max(pFloor, pB); // Safe Guard
	double cB = std::sqrt(gamma*pB/rhoP);

	double pP = gm1*(rhoEP - 0.5*rhoP*(uP*uP + vP*vP));
	pP = std::max(pFloor, pP); // Safe Guard
	double cP = std::sqrt(gamma*pP/rhoP);

	double uNP = uP*n(0) + vP*n(1);
	double JP = uNP + 2*cP/gm1;
	double uNB = JP - 2*cB/gm1;
	Eigen::Matrix<double, 4, 1> F;
	F(0) = 0.0;
	// F(1) = pB*n(0);
	// F(2) = pB*n(1);
	F(1) = pP*n(0);
	F(2) = pP*n(1);
	F(3) = 0.0;
        return F;
    }
};
