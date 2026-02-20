#pragma once
#include "boundaryFlux.hpp"
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <limits>

class inletFlux : public boundaryFlux {
private:
    double rho0_;
    double a0_;
    double alpha_;
    double t_;
    bool   transient_;

public:
    inletFlux(double rho0, double a0, double alpha, double t, bool transient)
        : rho0_(rho0), a0_(a0), alpha_(alpha), t_(t), transient_(transient) {}

    Eigen::Vector4d operator()(
        const Eigen::Vector4d& UP,
        double gamma,
        const Eigen::Vector2d& n
    ) const override {

	const double gm1 = gamma - 1.0;
	constexpr double rhoFloor = 1e-10;
	constexpr double pFloor = 1e-10;

	const double rhoP = std::max(rhoFloor, UP(0));
	const double uP = UP(1)/rhoP;
	const double vP = UP(2)/rhoP;
	const double rhoEP = UP(3);

	double pP = gm1*(rhoEP - 0.5*rhoP*(uP*uP + vP*vP));
	pP = std::max(pFloor, pP);

	const double cP = std::sqrt(gamma*pP/rhoP);
	Eigen::Matrix<double, 4, 1> F;
	const auto fillFlux = [&](double rho, double u, double v, double p, double rhoE) {
		F(0) = rho*(u*n(0) + v*n(1));
		F(1) = (rho*u*u + p)*n(0) + rho*u*v*n(1);
		F(2) = rho*u*v*n(0) + (rho*v*v + p)*n(1);
		F(3) = u*(rhoE + p)*n(0) + v*(rhoE + p)*n(1);
	};
	const double uNP = uP*n(0) + vP*n(1);

	// Robustness: if flow leaves domain through inlet boundary, extrapolate interior state.
	if (uNP > 0.0) {
		fillFlux(rhoP, uP, vP, pP, rhoEP);
		return F;
	}

	double pT = rho0_*a0_*a0_/gamma;
	double RTT = pT/rho0_;
	double JP = uNP + 2*cP/gm1;
	double dn = n(0)*std::cos(alpha_) + n(1)*std::sin(alpha_);
	
	// Safe Guard
	double A = gamma*RTT*dn*dn - 0.5*gm1*JP*JP;
	if (std::abs(A) < 1e-14) {
		// fallback: pick a small Mach or use interior Mach estimate
		// simplest:
		A = (A >= 0 ? 1e-14 : -1e-14);
	}

	double B = 4*gamma*RTT*dn/gm1;
	double C = 4*gamma*RTT/(gm1*gm1) - JP*JP;

	// Safe Guard
	double disc = B*B - 4*A*C;
	disc = std::max(0.0, disc);
	
	double MB = 0.1;
	if (std::abs(A) > 1e-14) {
		const double sqrtDisc = std::sqrt(disc);
		const double MB1 = (-B - sqrtDisc) / (2*A);
		const double MB2 = (-B + sqrtDisc) / (2*A);
		const bool m1Pos = MB1 > 0.0;
		const bool m2Pos = MB2 > 0.0;
		if (m1Pos && m2Pos) {
			MB = std::min(MB1, MB2);
		} else if (m1Pos) {
			MB = MB1;
		} else if (m2Pos) {
			MB = MB2;
		}
	}

	double RTB = RTT/(1 + 0.5*gm1*MB*MB);
	double pB = pT*std::pow(RTB/RTT, gamma/gm1);
	// double rhoB = pB/RTB;
	double rhoB = std::max(rhoFloor, pB/RTB);	
	double cB = std::sqrt(gamma*pB/rhoB);
	// double uNB = JP - 2*cB/gm1;
	double uB = MB*cB*std::cos(alpha_);
	double vB = MB*cB*std::sin(alpha_);
	double rhoEB = pB/gm1 + 0.5*rhoB*(uB*uB + vB*vB);

	// Eigen::Matrix<double, 4, 1> F;
	// F(0) = rhoB*(uB*n(0) + vB*n(1));
	// F(1) = (rhoB*uB*uB + pB)*n(0) + rhoB*uB*vB*n(1);
	// F(2) = rhoB*uB*vB*n(0) + (rhoB*vB*vB + pB)*n(1);
	// F(3) = uB*(rhoEB + pB)*n(0) + vB*(rhoEB + pB)*n(1);
	fillFlux(rhoB, uB, vB, pB, rhoEB);
	
        return F;
    }
};
