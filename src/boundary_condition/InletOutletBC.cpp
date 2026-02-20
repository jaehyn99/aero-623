#include "InletOutletBC.h"
#include "InletBC.h"
#include "OutletBC.h"

InletOutletBC::InletOutletBC(double rho0, double a0, double alpha, double pB, double gamma, bool transient):
    _inlet(std::make_unique<InletBC>(rho0, a0, alpha, gamma, transient)),
    _outlet(std::make_unique<OutletBC>(pB, gamma))
{}

Eigen::Vector4d InletOutletBC::computeFlux(const Eigen::Vector4d& UP, const Eigen::Vector2d& n) const{
    if (UP[1]*n[0] + UP[2]*n[1] < 0) return _inlet->computeFlux(UP, n); // Inflow
    return _outlet->computeFlux(UP, n); // Outflow
}

Eigen::Vector4d InletOutletBC::computeBoundaryState(const Eigen::Vector4d& UP, const Eigen::Vector2d& n) const{
    if (UP[1]*n[0] + UP[2]*n[1] < 0) return _inlet->computeBoundaryState(UP, n); // Inflow
    return _outlet->computeBoundaryState(UP, n); // Outflow
}

bool InletOutletBC::isTransient() const noexcept{ return _inlet->isTransient(); }
void InletOutletBC::setTransient(bool transient){ _inlet->setTransient(transient); }
void InletOutletBC::setTransientTime(double t) const{ _inlet->setTransientTime(t); }
void InletOutletBC::setTransientRho(double y) const{ _inlet->setTransientRho(y); }
void InletOutletBC::reset() const noexcept{ _inlet->reset(); }