#ifndef SECANT_SOLVER_H
#define SECANT_SOLVER_H

#include "NonlinearSolver.h"

class SecantSolver : public NonlinearSolver<1, 1, double>{
    public:
    using typename NonlinearSolver<1, 1, double>::DomainType;
    using typename NonlinearSolver<1, 1, double>::RangeType;
    using typename NonlinearSolver<1, 1, double>::Function;

    explicit SecantSolver(const Function& func, double ftol=1.0e-6, double xtol=1.0e-6, std::size_t maxIter=30):
        NonlinearSolver<1, 1, double>(func, ftol, xtol, maxIter) {};
    
    NLStatus solve(DomainType& x0, double&& x1) const noexcept override{
        double f0 = this->_f(x0);
        double f1 = this->_f(x1);

        for (std::size_t iter = 0; iter < this->_maxIter; iter++){
            if (std::abs(f0-f1) < _ftol) return NLStatus::SingularityError;

            double tempx = x0;
            double tempf = f0;
            double dx = f0 * (x0-x1)/(f0-f1);
            x0 -= dx;
            f0 = this->_f(x0);
            if (this->inputConverged(x0, dx) || this->outputConverged(f0)) return NLStatus::Success;

            x1 = tempx;
            f1 = tempf;
        }
        return NLStatus::NoConvergence;
    }
};

#endif