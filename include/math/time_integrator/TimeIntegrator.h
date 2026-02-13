#ifndef TIME_INTEGRATOR_H
#define TIME_INTEGRATOR_H

#include <functional>
#include "StateMesh.h"

class TimeIntegrator{
    public:
    using Function = std::function<StateMesh(double, const StateMesh&)>;
    virtual ~TimeIntegrator() = default;
    virtual StateMesh integrate(const Function& f, const StateMesh& u0, double t, double dt) const = 0; // Global time-stepping
    virtual StateMesh integrate(const Function& f, const StateMesh& u0, double t, const Eigen::ArrayXd& dt) const = 0; // Local time-stepping
};

#endif