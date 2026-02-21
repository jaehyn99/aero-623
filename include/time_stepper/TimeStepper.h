#ifndef TIME_STEPPER_H
#define TIME_STEPPER_H

#include "Eigen/Dense"

class StateMesh;
class FVFlux;
class TimeStepper{
    public:
    TimeStepper(double minCFL, double gamma, std::shared_ptr<FVFlux> flux): _minCFL(minCFL), _gamma(gamma), _flux(flux) {}
    virtual ~TimeStepper() = default;
    // u: state vector, s: wave-speed vector on the edges
    virtual Eigen::ArrayXd dt(const StateMesh& u) const noexcept = 0;

    protected:
    double _minCFL;
    double _gamma;
    std::shared_ptr<FVFlux> _flux;
};

#endif