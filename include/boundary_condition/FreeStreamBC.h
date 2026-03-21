#ifndef FREE_STREAM_BC_H
#define FREE_STREAM_BC_H

#include "BoundaryCondition.h"

class FreeStreamBC : public BoundaryCondition{
	public:
	FreeStreamBC(double gamma): _gamma(gamma) {}
	Eigen::Vector4d computeFlux(const Eigen::Vector4d& UP, const Eigen::Vector2d& n) const override;
	Eigen::Vector4d computeBoundaryState(const Eigen::Vector4d& UP, const Eigen::Vector2d& n) const override{ return UP; };

	protected:
	double _gamma;
};

#endif