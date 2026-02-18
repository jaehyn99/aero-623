#ifndef INLET_H
#define INLET_H

#include "BoundaryCondition.h"

class InletBC: public BoundaryCondition {
	public:
    InletBC(double rho0, double a0, double alpha, double gamma/*, double t, bool transient*/);
	Eigen::Vector4d computeFlux(const Eigen::Vector4d& UP, const Eigen::Vector2d& n) const override;
    Eigen::Vector4d computeBoundaryState(const Eigen::Vector4d& UP, const Eigen::Vector2d& n) const override;

	protected:
    double _rho0;
    double _a0;
    double _alpha;
	double _gamma;
    // double _t;
    // bool   _transient;
};

#endif