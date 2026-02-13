#ifndef GLOBAL_TIME_STEPPER_H
#define GLOBAL_TIME_STEPPER_H

#include "TimeStepper.h"

class GlobalTimeStepper: public TimeStepper{
    public:
    GlobalTimeStepper(double minCFL): TimeStepper(minCFL) {}
    Eigen::ArrayXd dt(const StateMesh& u, const Eigen::ArrayXd& s) const noexcept override;
}

#endif