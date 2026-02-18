#ifndef FV_SOLVER_H
#define FV_SOLVER_H

#include "Eigen/Dense"

class FVResidual;
class StateMesh;
class TimeIntegrator;
class TimeStepper;
class FVSolver{
    public:
    FVSolver(std::shared_ptr<FVResidual>, std::shared_ptr<TimeIntegrator>, std::shared_ptr<TimeStepper>);
    FVSolver(FVSolver&&);
    FVSolver& operator=(FVSolver&&);
    ~FVSolver();

    virtual void solve(StateMesh&) const = 0;
    auto getResult() const noexcept{ return std::move(_result); }
    auto getNorm() const noexcept{ return std::move(_l1norm); }

    protected:
    std::shared_ptr<FVResidual> _residual;
    std::shared_ptr<TimeIntegrator> _integrator;
    std::shared_ptr<TimeStepper> _stepper;
    mutable std::vector<Eigen::MatrixXd> _result;
    mutable std::vector<double> _l1norm;
};

#endif