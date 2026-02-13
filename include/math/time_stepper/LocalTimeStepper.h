#ifndef LOCAL_TIME_STEPPER_H
#define LOCAL_TIME_STEPPER_H

#include "TimeStepper.h"

class LocalTimeStepper: public TimeStepper{
    public:
    LocalTimeStepper(double minCFL): TimeStepper(minCFL) {}
    Eigen::ArrayXd dt(const StateMesh& u, const Eigen::ArrayXd& s) const noexcept override;
}

#endif