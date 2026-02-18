#ifndef TIME_INTEGRATOR_H
#define TIME_INTEGRATOR_H

#include <functional>
#include "Eigen/Dense"

class StateMesh;
class TimeIntegrator{
    public:
    using Function = std::function<Eigen::MatrixXd(double, Eigen::MatrixXd&)>;
    virtual ~TimeIntegrator() = default;
    virtual void integrate(const Function& f, StateMesh& u0, double t, double dt) const = 0; // Global time-stepping
    virtual void integrate(const Function& f, StateMesh& u0, double t, const Eigen::ArrayXd& dt) const = 0; // Local time-stepping
};

#endif