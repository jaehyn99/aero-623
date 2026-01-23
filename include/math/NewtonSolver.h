#ifndef NEWTON_SOLVER_H
#define NEWTON_SOLVER_H

#include "NonlinearSolver.h"
#include "Derivative.h"

template<int N, int M = N>
class NewtonSolver : public NonlinearSolver<N, M>{
    public:
    using typename NonlinearSolver<N, M>::DomainType;
    using typename NonlinearSolver<N, M>::RangeType;
    using typename NonlinearSolver<N, M>::Function;
    using DerivativeType = std::conditional_t<M == 1, double, Eigen::Matrix<double, M, N>>;
    using DFunction = std::function<DerivativeType(const DomainType&)>;

    explicit NewtonSolver(const Function& func, const DFunction& dfunc=nullptr,
                          double ftol=1.0e-6, double xtol=1.0e-6, std::size_t maxIter=20):
        NonlinearSolver<N, M>(func, ftol, xtol, maxIter),
        _df(dfunc)
    {
        if (!_df){ // Derivative not provided, use numerical differentiation
            if constexpr(N == 1 && M == 1) _df = [this](const DomainType& x){ return df(this->_f, x); };
            else _df = [this](const DomainType& x){ return Jacobian(this->_f, x); };
        }
    }

    NLStatus solve(DomainType& x) const noexcept override{
        if constexpr(N == 1){
            for (std::size_t iter = 0; iter < this->_maxIter; iter++){
                DerivativeType dfx = _df(x);
                if (dfx == 0.0) return NLStatus::SingularityError;
                DomainType dx = this->_f(x)/dfx;
                x -= dx;
                RangeType fx = this->_f(x);
                if (this->inputConverged(x, dx) || this->outputConverged(fx)) return NLStatus::Success;
            }
            return NLStatus::NoConvergence;
        }

        else{
            RangeType fx = this->_f(x);
            DerivativeType Jx = _df(x);
            if (N == Eigen::Dynamic || M == Eigen::Dynamic){
                if (fx.size() < x.size() || fx.size() < Jx.rows() || fx.size() != Jx.cols()) return NLStatus::InvalidArgument;
            }
            for (std::size_t iter = 0; iter < this->_maxIter; iter++){
                Eigen::ColPivHouseholderQR<Eigen::Ref<DerivativeType>> qr(Jx);
                if (qr.rank() < fx.size()) return NLStatus::SingularityError;
                auto dx = qr.solve(fx);
                x -= dx;
                fx = this->_f(x);
                if (this->inputConverged(x, dx) || this->outputConverged(fx)) return NLStatus::Success;
                Jx = _df(x);
            }
            return NLStatus::NoConvergence;
        }
    }
    
    void setFunction(const Function& f) const noexcept override{
        this->_f = f;
        if constexpr (N == 1) _df = [this](const DomainType& x){ return df(this->_f, x); };
        else _df = [this](const DomainType& x){ return Jacobian(this->_f, x); };
    }

    void setFunction(const Function& f, const DFunction& df) const noexcept{
        this->_f = f;
        _df = df;
    }

    protected:
    mutable DFunction _df;
};

#endif