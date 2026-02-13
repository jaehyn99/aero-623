#ifndef FV_SOLVER_H
#define FV_SOLVER_H

#include "Eigen/Dense"

class Residual;
class StateMesh;
class TimeIntegrator;
class FVSolver{
    public:
    FVSolver(std::unique_ptr<Residual>, std::unique_ptr<TimeIntegrator>, bool, double=1e-5, double=1000);
    FVSolver(FVSolver&&);
    FVSolver& operator=(FVSolver&&);
    ~FVSolver();

    void solve(StateMesh&) const noexcept;
    auto getResult() const noexcept{ return std::move(_result); }
    auto getNorm() const noexcept{ return std::move(_l1norm); }

    protected:
    std::unique_ptr<Residual> _residual;
    std::unique_ptr<TimeIntegrator> _integrator;
    bool _isSteady;
    double _tol; // for steady solve
    double _tMax; // for unsteady solve
    mutable std::vector<Eigen::ArrayXd> _result;
    mutable std::vector<double> _l1norm;
};

#endif