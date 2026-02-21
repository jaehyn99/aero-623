#include "OutletBC.h"

Eigen::Vector4d OutletBC::computeFlux(const Eigen::Vector4d& UP, const Eigen::Vector2d& n) const{
	double gm1 = _gamma - 1;
    
    double rhoP = UP(0); 
	double uP = UP(1)/rhoP;
	double vP = UP(2)/rhoP;
	double rhoEP = UP(3);
	
    double pP = gm1*(rhoEP - 0.5*rhoP*(uP*uP + vP*vP));
	double SP = pP/std::pow(rhoP, _gamma);
	double rhoB = std::pow(_pB/SP, 1/_gamma);
	double cB = std::sqrt(_gamma*_pB/rhoB);
	double uNP = uP*n(0) + vP*n(1);
	double cP = std::sqrt(_gamma*pP/rhoP);
	double JP = uNP + 2*cP/gm1;
	double uNB = JP - 2*cB/gm1; // Boundary normal velocity

	double uB = uP - n(0)*(uP*n(0) + vP*n(1) - uNB);
	double vB = vP - n(1)*(uP*n(0) + vP*n(1) - uNB);
	// Supersonic outflow, use interior states
	if (uB*uB + vB*vB > cB*cB) return {rhoP*(uP*n(0) + vP*n(1)),
									   (rhoP*uP*uP + pP)*n(0) + rhoP*uP*vP*n(1),
									   rhoP*uP*vP*n(0) + (rhoP*vP*vP + pP)*n(1),
									   uP*(rhoEP + pP)*n(0) + vP*(rhoEP + pP)*n(1)};

	// Subsonic outflow, use boundary states
	double rhoEB = _pB/gm1 + 0.5*rhoB*(uB*uB + vB*vB);
	return {rhoB*(uB*n(0) + vB*n(1)),
			(rhoB*uB*uB + _pB)*n(0) + rhoB*uB*vB*n(1),
			rhoB*uB*vB*n(0) + (rhoB*vB*vB + _pB)*n(1),
			uB*(rhoEB + _pB)*n(0) + vB*(rhoEB + _pB)*n(1)};
}

Eigen::Vector4d OutletBC::computeBoundaryState(const Eigen::Vector4d& UP, const Eigen::Vector2d& n) const{
	double gm1 = _gamma - 1;
    
    double rhoP = UP(0); 
	double uP = UP(1)/rhoP;
	double vP = UP(2)/rhoP;
	double rhoEP = UP(3);
	
    double pP = gm1*(rhoEP - 0.5*rhoP*(uP*uP + vP*vP));
	double SP = _pB/std::pow(rhoP, _gamma);
	double rhoB = std::pow(_pB/SP, 1/_gamma);
	double cB = std::sqrt(_gamma*_pB/rhoB);

	double uNP = uP*n(0) + vP*n(1);
	double cP = std::sqrt(_gamma*pP/rhoP);
	double JP = uNP + 2*cP/gm1;
	double uNB = JP - 2*cB/gm1;
	double uB = uP - n(0)*(uP*n(0) + vP*n(1) - uNB);
	double vB = vP - n(1)*(uP*n(0) + vP*n(1) - uNB);
	return {rhoB, rhoB*uB, rhoB*vB, cB};
}