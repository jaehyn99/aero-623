#ifndef OUTLET_H
#define OUTLET_H

#include "BoundaryCondition.h"

class OutletBC : public BoundaryCondition {
	public:
    OutletBC(double pB, double gamma): _pB(pB), _gamma(gamma) {}
	Eigen::Vector4d computeFlux(const Eigen::Vector4d& UP, const Eigen::Vector2d& n) const override;
	Eigen::Vector4d computeBoundaryState(const Eigen::Vector4d& UP, const Eigen::Vector2d& n) const override;

	protected:
	double _pB;
	double _gamma;
};

#endif