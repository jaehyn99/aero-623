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
        Eigen::ArrayXd dt = _stepper->dt(u);
        _integrator->integrate(func, u.matrix(), 0, dt);

        norm = func(0, u.matrix()).lpNorm<1>();      
        _l1norm.push_back(norm);
        _result.emplace_back(u.matrix());
        isConverged = norm/_l1norm.front() <= _tol;
        // std::cout << norm << std::endl;

        double maxP = 0;
        for (Eigen::Index e = 0; e < u.cellCount(); e++){
            double rhoE = u(3, e);
            double KE = 0.5*(u(1,e)*u(1,e) + u(2,e)*u(2,e)) / u(0,e);
            double P = 0.4 * (rhoE - KE);
            if (P > maxP) maxP = P;
        }
    }
}