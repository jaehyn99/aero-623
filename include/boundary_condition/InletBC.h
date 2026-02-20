#ifndef INLET_H
#define INLET_H

#include "BoundaryCondition.h"

class InletBC: public BoundaryCondition {
	public:
    InletBC(double rho0, double a0, double alpha, double gamma, bool transient=false);
	Eigen::Vector4d computeFlux(const Eigen::Vector4d& UP, const Eigen::Vector2d& n) const override;
    Eigen::Vector4d computeBoundaryState(const Eigen::Vector4d& UP, const Eigen::Vector2d& n) const override;

    bool isTransient() const noexcept{ return _transient; }
    void setTransient(bool transient){ _transient = transient; }
    void setTransientTime(double t) const;
    void setTransientRho(double y) const;
    void reset() const noexcept; // sets back to steady-state value

	protected:
    double _rho0;
    mutable double _rhoTrans; // transient value of _rho0;
    double _a0;
    double _alpha;
	double _gamma;
    mutable double _t;
    bool _transient;
};

#endif