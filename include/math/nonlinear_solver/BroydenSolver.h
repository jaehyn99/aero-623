#ifndef BROYDEN_SOLVER_H
#define BROYDEN_SOLVER_H

#include "NonlinearSolver.h"

template<int N>
using BroydenSolverBase = NonlinearSolver<N, N, Eigen::Matrix<double, N, N>&>;

template<int N>
class BroydenSolver : public BroydenSolverBase<N>{
    static_assert(N == Eigen::Dynamic || N > 1, "For multi-dimensional problems only. Use secant method for 1-D problems.");
    public:
    using typename BroydenSolverBase<N>::DomainType;
    using typename BroydenSolverBase<N>::Function;
    using DerivativeType = Eigen::Matrix<double, N, N>&;

    explicit BroydenSolver(const Function& func, bool good=true, double ftol=1.0e-6, double xtol=1.0e-6, std::size_t maxIter=20):
        BroydenSolverBase<N>(func, ftol, xtol, maxIter),
        _good(good)
    {}

    NLStatus solve(DomainType& x, DerivativeType& Jinv) const noexcept override{
        DomainType f = this->_f(x);
        if (N == Eigen::Dynamic){
            if (x.size() == 1 || f.size() != x.size () || f.size() != Jinv.rows() || f.size() != Jinv.cols())
                return NLStatus::InvalidArgument;
        }

        for (std::size_t iter = 0; iter < this->_maxIter; iter++){
            DomainType dx = Jinv*f;
            x -= dx;
            DomainType f1 = this->_f(x);
            DomainType df = f1 - f;
            if (this->inputConverged(x, dx) || this->outputConverged(f1)) return NLStatus::Success;
            
            // Update the inverse Jacobian and the function value
            if (_good){ // "Good" Broyden method
                double denom = dx.dot(Jinv*df);
                if (denom == 0.0) return NLStatus::SingularityError;
                Jinv += ((dx - Jinv*df).outer(dx) / denom) * Jinv;
            } else{ // "Bad" Broyden method
                double denom = df.squaredNorm();
                if (denom == 0.0) return NLStatus::SingularityError;
                Jinv += (dx - Jinv*df).outer(df) / denom;
            }
            f = f1;
        }
        return NLStatus::NoConvergence;
    }

    protected:
    bool _good; // whether to use "good" or "bad" Broyden method
};

#endif