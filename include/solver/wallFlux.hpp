#pragma once
#include "boundaryFlux.hpp"
#include <Eigen/Dense>
#include <algorithm>

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
        const double gm1 = gamma - 1.0;

        // Minimal robustness guards
        constexpr double rhoFloor = 1e-10;
        constexpr double pFloor = 1e-10;

        const double rhoP = std::max(rhoFloor, UP(0));
        const double uP = UP(1) / rhoP;
        const double vP = UP(2) / rhoP;
        const double rhoEP = UP(3);

        // Slip-wall tangential velocity state (remove normal velocity)
        const double un = uP * n(0) + vP * n(1);
        const double uB = uP - un * n(0);
        const double vB = vP - un * n(1);

        // Wall pressure from slip-wall projected velocity
        double pB = gm1 * (rhoEP - 0.5 * rhoP * (uB * uB + vB * vB));
        pB = std::max(pFloor, pB);

        Eigen::Vector4d F;
        F(0) = 0.0;
        F(1) = pB * n(0);
        F(2) = pB * n(1);
        F(3) = 0.0;
        return F;
    }
};
