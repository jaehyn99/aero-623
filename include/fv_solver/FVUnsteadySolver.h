#ifndef FV_UNSTEADY_SOLVER_H
#define FV_UNSTEADY_SOLVER_H

#include "Solver.h"
class FVUnSteadySolver: public Solver{
    public:
    FVUnSteadySolver(std::shared_ptr<Residual>, std::shared_ptr<TimeIntegrator>, std::shared_ptr<TimeStepper>, std::size_t, std::size_t);
    void solve(StateMesh&) const override;

    protected:
    std::size_t _saveEveryNIterations;
    std::size_t _maxIterations;
};

#endif