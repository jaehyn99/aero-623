#include "FVUnsteadySolver.h"
#include "FVResidual.h"
#include "InletBC.h"
#include "InletOutletBC.h"
#include "LocalTimeStepper.h"
#include "StateMesh.h"
#include "TimeIntegrator.h"
#include "TimeStepper.h"
#include <iostream>

FVUnSteadySolver::FVUnSteadySolver(std::shared_ptr<FVResidual> residual, std::shared_ptr<TimeIntegrator> integrator,
                               std::shared_ptr<TimeStepper> stepper, std::size_t saveEveryIterations, std::size_t maxIterations):
    FVSolver(residual, integrator, stepper),
    _saveEveryIterations(saveEveryIterations),
    _maxIterations(maxIterations)
{}

void FVUnSteadySolver::solve(StateMesh& u) const{
    auto func = [this, &u](double t, Eigen::MatrixXd& x)
    {
        for (Eigen::Index i = 0; i < u.bcCount(); i++){
            BoundaryCondition* bcPtr = u.bc(i).get();
            if (auto inlet = dynamic_cast<InletBC*>(bcPtr)){
                if (inlet->isTransient()) inlet->setTransientTime(t);
            }
            else if (auto inlet = dynamic_cast<InletOutletBC*>(bcPtr)){
                if (inlet->isTransient()) inlet->setTransientTime(t);
            }
        }
        u.matrix() = std::move(x);
        return _residual->computeResidual(u);
    };
    double norm = func(0, u.matrix()).lpNorm<1>();
    _l1norm.push_back(norm); // the first L1-norm

    _result.reserve(_maxIterations/_saveEveryIterations+1);
    for (std::size_t iter = 0; iter < _maxIterations; iter++){
        double dt = _stepper->dt(u).minCoeff();
        _integrator->integrate(func, u.matrix(), 0, dt);

        norm = func(0, u.matrix()).lpNorm<1>();      
        _l1norm.push_back(norm);
        if (iter % _saveEveryIterations == 0) _result.emplace_back(u.matrix());
    }
    if (_maxIterations % _saveEveryIterations != 0) _result.emplace_back(u.matrix());
}