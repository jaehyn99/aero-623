// #include "FVSolver.h"
// #include "GlobalTimeStepper.h"
// #include "LocalTimeStepper.h"
// #include "Residual.h"
// #include "StateMesh.h"
// #include "TimeIntegrator.h"

// FVSolver::FVSolver(std::unique_ptr<Residual> residual, std::unique_ptr<TimeIntegrator> integrator, bool isSteady, double tol, double tMax):
//     _residual(std::move(residual)),
//     _integrator(std::move(integrator)),
//     _isSteady(isSteady),
//     _tMax(tMax)
// {
//     assert(_residual && _integrator && _tMax > 0);
// }

// FVSolver::FVSolver(FVSolver&&) = default;
// FVSolver& operator=(FVSolver&&) = default; 
// FVSolver::~FVSolver() = default;

// void FVSolver::solve(StateMesh& u) const noexcept{
//     auto func = [auto F = std::move(_residual)](double, const StateMesh& u) -> Eigen::ArrayXd
//         { return F->calculateFlux(u); };
//     _l1norm.push_back(func(0, u).lpNorm<1>); // the norm to be reached in order to consider the solve converged

//     if (_isSteady){
//         LocalTimeStepper stepper(1.0);
//         Eigen::ArrayXd s; // where does this come from?
//         auto dt = stepper.dt(u, s);

//         bool isConverged = false;
//         while (!isConverged){
//             u = _integrator->integrate(func, u, 0, dt);
//             _l1norm.push_back(func(0, u).lpNorm<1>);
//             isConverged = _l1norm.back()/_l1norm.front() <= _tol;
//         }
//         _result.emplace_back(std::move(u));
//     } else{
//         GlobalTimeStepper stepper(1.0);
//         Eigen::ArrayXd s; // where does this come from?
//         double dt = stepper.dt(u, s)[0]; // just a single time step for all elements
//         double t = 0;
//         _result.reserve(static_cast<int>(_tMax/dt));
//         _result.push_back(u);
//         _l1norm.reserve(static_cast<int>(_tMax/dt));

//         while (t < _tMax){
//             _result.push_back(_integrator->integrate(func, u, t, dt));
//             _l1norm.push_back(func(t, _result.back()).lpNorm<1>);
//             t += dt;
//         }
//     }
// }