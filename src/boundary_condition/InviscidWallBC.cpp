#include "InviscidWallBC.h"

Eigen::Vector4d InviscidWallBC::computeFlux(const Eigen::Vector4d& UP, const Eigen::Vector2d& n) const{
	double gm1 = _gamma - 1;

	double rhoP = UP(0);
	double uP = UP(1)/rhoP;
	double vP = UP(2)/rhoP;
	double rhoEP = UP(3);

	double uB = uP - n(0)*(uP*n(0) + vP*n(1));
	double vB = vP - n(1)*(uP*n(0) + vP*n(1));
	double pB = gm1*(rhoEP - 0.5*rhoP*(uB*uB + vB*vB));
	return {0.0, pB*n(0), pB*n(1), 0.0};
}

Eigen::Vector4d InviscidWallBC::computeBoundaryState(const Eigen::Vector4d& UP, const Eigen::Vector2d& n) const{
	double gm1 = _gamma - 1;

	double rhoP = UP(0);
	double uP = UP(1)/rhoP;
	double vP = UP(2)/rhoP;
	double rhoEP = UP(3);

	double uB = uP - n(0)*(uP*n(0) + vP*n(1));
	double vB = vP - n(1)*(uP*n(0) + vP*n(1));
	double pB = gm1*(rhoEP - 0.5*rhoP*(uB*uB + vB*vB));	
	double cB = std::sqrt(_gamma*pB/rhoP);
	return {rhoP, rhoP*uB, rhoP*vB, cB};
}