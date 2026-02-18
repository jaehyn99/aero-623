#ifndef WALL_H
#define WALL_H

#include "BoundaryCondition.h"

class InviscidWallBC: public BoundaryCondition{
	public:
	InviscidWallBC(double gamma): _gamma(gamma) {}
	Eigen::Vector4d computeFlux(const Eigen::Vector4d& UP, const Eigen::Vector2d& n) const override;
	Eigen::Vector4d computeBoundaryState(const Eigen::Vector4d& UP, const Eigen::Vector2d& n) const override;

	protected:
	double _gamma;
};

#endif