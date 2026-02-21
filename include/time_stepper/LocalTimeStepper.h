#ifndef LOCAL_TIME_STEPPER_H
#define LOCAL_TIME_STEPPER_H

#include "TimeStepper.h"

class LocalTimeStepper: public TimeStepper{
    public:
    LocalTimeStepper(double minCFL, double gamma, std::shared_ptr<FVFlux> flux): TimeStepper(minCFL, gamma, flux) {}
    Eigen::ArrayXd dt(const StateMesh& u) const noexcept override;
};

#endif