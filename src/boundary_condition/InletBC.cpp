#include "InletBC.h"
#include <iostream>

InletBC::InletBC(double rho0, double a0, double alpha, double gamma/*, double t, bool transient*/):
    _rho0(rho0),
    _a0(a0),
    _alpha(alpha),
    _gamma(gamma)
    // _t(t),
    // _transient(transient)
{}

Eigen::Vector4d InletBC::computeFlux(const Eigen::Vector4d& UP, const Eigen::Vector2d& n) const{
	double gm1 = _gamma - 1;
    
    double rhoP = UP(0);
	double uP = UP(1)/rhoP;
	double vP = UP(2)/rhoP;
	double rhoEP = UP(3);
	
    double pP = gm1*(rhoEP - 0.5*rhoP*(uP*uP + vP*vP));
	double cP = std::sqrt(_gamma*pP/rhoP);
	double uNP = uP*n(0) + vP*n(1);
	double pT = _rho0*_a0*_a0/_gamma;
	double RTT = pT/_rho0;
	double JP = uNP + 2*cP/gm1;
	double dn = n(0)*std::cos(_alpha) + n(1)*std::sin(_alpha);

	double A = _gamma*RTT*dn*dn - 0.5*gm1*JP*JP;
	double B = 4*_gamma*RTT*dn/gm1;
	double C = 4*_gamma*RTT/(gm1*gm1) - JP*JP;
	double MB;
	if (B*B - 4*A*C >= 0){
		double MB1 = (-B + std::sqrt(B*B - 4*A*C))/(2*A);
		double MB2 = (-B - std::sqrt(B*B - 4*A*C))/(2*A);
		// if (MB1*MB2 < 0) MB = std::max(MB1, MB2);
		// else if (MB1 > 0) MB = std::min(MB1, MB2);
		// else throw std::runtime_error("ERROR: Both Mach numbers are negative.");
		MB = MB1*MB2 < 0 ? std::max(MB1, MB2) : std::min(std::abs(MB1), std::abs(MB2));
	} else MB = std::sqrt(uP*uP + vP*vP) / cP; // take the interior Mach number value
	// } else throw std::runtime_error("ERROR: Negative discriminant. Cannot solve for the Mach number.");

	double RTB = RTT/(1 + 0.5*gm1*MB*MB);
	double pB = pT*std::pow(RTB/RTT, _gamma/gm1);
	double rhoB = pB/RTB;
	double cB = std::sqrt(_gamma*pB/rhoB);
	double uB = MB*cB*std::cos(_alpha);
	double vB = MB*cB*std::sin(_alpha);
	double rhoEB = pB/gm1 + 0.5*rhoB*(uB*uB + vB*vB);

	return {rhoB*(uB*n(0) + vB*n(1)),
            (rhoB*uB*uB + pB)*n(0) + rhoB*uB*vB*n(1),
            rhoB*uB*vB*n(0) + (rhoB*vB*vB + pB)*n(1),
            uB*(rhoEB + pB)*n(0) + vB*(rhoEB + pB)*n(1)};
}

Eigen::Vector4d InletBC::computeBoundaryState(const Eigen::Vector4d& UP, const Eigen::Vector2d& n) const{
	double gm1 = _gamma - 1;
    double rhoP = UP(0); // Will make this a function based on the transient flag and time soon
	double uP = UP(1)/rhoP;
	double vP = UP(2)/rhoP;
	double rhoEP = UP(3);
	
    double pP = gm1*(rhoEP - 0.5*rhoP*(uP*uP + vP*vP));
	double cP = std::sqrt(_gamma*pP/rhoP);
	double uNP = uP*n(0) + vP*n(1);
	double pT = _rho0*_a0*_a0/_gamma;
	double RTT = pT/_rho0;
	double JP = uNP + 2*cP/gm1;
	double dn = n(0)*std::cos(_alpha) + n(1)*std::sin(_alpha);

	double A = _gamma*RTT*dn*dn - 0.5*gm1*JP*JP;
	double B = 4*_gamma*RTT*dn/gm1;
	double C = 4*_gamma*RTT/(gm1*gm1) - JP*JP;
	double MB1 = (-B + std::sqrt(B*B - 4*A*C))/(2*A);
	double MB2 = (-B - std::sqrt(B*B - 4*A*C))/(2*A);
	double MB = MB1*MB2 < 0 ? std::max(MB1, MB2) : std::min(std::abs(MB1), std::abs(MB2));

	double RTB = RTT/(1 + 0.5*gm1*MB*MB);
	double pB = pT*std::pow(RTB/RTT, _gamma/gm1);
	double rhoB = pB/RTB;
	double cB = std::sqrt(_gamma*pB/rhoB);
	double uB = MB*cB*std::cos(_alpha);
	double vB = MB*cB*std::sin(_alpha);
	return {rhoB, rhoB*uB, rhoB*vB, cB};
}