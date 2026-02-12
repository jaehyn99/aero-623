#ifndef TIME_STEPPER_H
#define TIME_STEPPER_H

#include "Eigen/Dense"

class StateMesh;
class TimeStepper{
    public:
    TimeStepper(double minCFL): _minCFL(minCFL) { assert(_minCFL > 0); }
    virtual ~TimeStepper() = default;
    // u: state vector, s: wave-speed vector on the edges
    virtual Eigen::ArrayXd dt(const StateMesh& u, const Eigen::ArrayXd& s) const noexcept = 0;

    protected:
    double _minCFL;
};

#endif