#ifndef FV_UNSTEADY_SOLVER_H
#define FV_UNSTEADY_SOLVER_H

#include "FVSolver.h"
class FVUnSteadySolver: public FVSolver{
    public:
    FVUnSteadySolver(std::shared_ptr<FVResidual>, std::shared_ptr<TimeIntegrator>, std::shared_ptr<TimeStepper>, std::size_t, std::size_t);
    void solve(StateMesh&) const override;

    protected:
    std::size_t _saveEveryIterations;
    std::size_t _maxIterations;
};

#endif