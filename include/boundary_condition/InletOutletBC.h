#ifndef INLET_OUTLET_BC_H
#define INLET_OUTLET_BC_H

#include "BoundaryCondition.h"
class InletBC;
class OutletBC;
class InletOutletBC: public BoundaryCondition {
	public:
    InletOutletBC(double rho0, double a0, double alpha, double pB, double gamma, bool transient=false);
    InletOutletBC(InletOutletBC&&);
    InletOutletBC& operator=(InletOutletBC&&);
    ~InletOutletBC();

	Eigen::Vector4d computeFlux(const Eigen::Vector4d& UP, const Eigen::Vector2d& n) const override;
    Eigen::Vector4d computeBoundaryState(const Eigen::Vector4d& UP, const Eigen::Vector2d& n) const override;

    bool isTransient() const noexcept;
    void setTransient(bool transient);
    void setTransientTime(double t) const;
    void setTransientRho(double y) const;
    void reset() const noexcept;

	protected:
    std::unique_ptr<InletBC> _inlet;
    std::unique_ptr<OutletBC> _outlet;
};

#endif