#include "FVSteadySolver.h"
#include "FVResidual.h"
#include "LocalTimeStepper.h"
#include "StateMesh.h"
#include "TimeIntegrator.h"
#include "TimeStepper.h"
#include <iostream>

FVSteadySolver::FVSteadySolver(std::shared_ptr<FVResidual> residual, std::shared_ptr<TimeIntegrator> integrator,
                               std::shared_ptr<TimeStepper> stepper, double tol):
    FVSolver(residual, integrator, stepper),
    _tol(tol)
{}

void FVSteadySolver::solve(StateMesh& u) const{
    auto func = [this, &u](double t, Eigen::MatrixXd& x)
    {
        u.matrix() = std::move(x);
        return _residual->computeResidual(u);
    };
    double norm = func(0, u.matrix()).lpNorm<1>();
    _l1norm.push_back(norm); // the first L1-norm

    bool isConverged = false;
    while (!isConverged){
        std::cout << norm << std::endl;
        double dt = _stepper->dt(u).minCoeff();
        _integrator->integrate(func, u, 0, dt);

        norm = func(0, u.matrix()).lpNorm<1>();      
        _l1norm.push_back(norm);
        isConverged = norm/_l1norm.front() <= _tol;
    }
    _result.emplace_back(u.matrix());
}