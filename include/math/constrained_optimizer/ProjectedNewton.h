#ifndef PROJECTED_NEWTON_H
#define PROJECTED_NEWTON_H

#include "ConstrainedOptimizer.h"
#include "Eigen/Cholesky"

template<int N, typename = std::enable_if_t<(N == Eigen::Dynamic || N >= 1)>>
class ProjectedNewton : public ConstrainedOptimizer<N>{
    public:
    using typename ConstrainedOptimizer<N>::DomainType;
    using typename ConstrainedOptimizer<N>::Function;
    using typename ConstrainedOptimizer<N>::GFunction;
    using HessianType = std::conditional_t<N == 1, double, Eigen::Matrix<double, N, N>>;
    using HFunction = std::function<HessianType(const DomainType&)>;

    explicit ProjectedNewton(const Function& func, const DomainType& lower, const DomainType& upper,
                             const GFunction& gfunc=nullptr, const HFunction& Hfunc=nullptr,
                             double ftol=1.0e-6, double gtol=1.0e-6, double xtol=1.0e-6, std::size_t maxIter=30):
        ConstrainedOptimizer<N>(func, gfunc, ftol, gtol, xtol, maxIter),
        _L(lower),
        _U(upper),
        _H(Hfunc)
    {
        if (!Hfunc){
            if constexpr(N == 1) _H = [this](const DomainType& x){ return df2(this->_f, x); };
            else _H = [this](const DomainType& x){ return Hessian(this->_f, x); };
        }
    }

    COStatus minimize(DomainType& x) const noexcept override{
        constexpr double c = 1e-4;
        double fx = this->_f(x);
        DomainType gx = this->_g(x);
        HessianType Hx = _H(x);

        if constexpr(N == 1){
            if (this->gradientConverged(gx) && Hx < 0.0) return COStatus::Success;
            for (std::size_t iter = 0; iter < this->_maxIter; iter++){
                // Clamp the step size and advance
                DomainType dx = (Hx > 0.0) ? -gx/Hx : -gx;
                double alpha = dx > 0 ? std::min(1.0, (_U-x)/dx) : std::min(1.0, (_L-x)/dx); // maximum step size before hitting a bound
                double fx1 = this->_f(x+alpha*dx);
                while (fx1 > fx + c*alpha*gx*dx){
                    alpha *= 0.5;
                    fx1 = this->_f(x+alpha*dx);
                }
                x += alpha*dx;
                gx = this->_g(x);
        
                // Convergence check
                if (this->inputConverged(x, dx) || this->gradientConverged(gx) || this->outputConverged(fx, fx1)) return COStatus::Success;
                fx = fx1;
                Hx = _H(x);
            }
            return COStatus::NoConvergence;
        }

        else{
            if constexpr (N == Eigen::Dynamic){
                if (gx.size() != Hx.rows() || gx.size() != Hx.cols()) return COStatus::InvalidArgument;
            }

            for (std::size_t iter = 0; iter < this->_maxIter; iter++){
                // Check if the Hessian is positive definite
                Eigen::LLT<Eigen::Ref<HessianType>> llt(Hx);
                double lambda = 1.0e-8; // Levenberg–Marquardt damping constant
                while (llt.info() != Eigen::Success){
                    Hx = llt.matrixL();
                    Hx = Hx * Hx.transpose();
                    for (Eigen::Index i = 0; i < gx.size(); i++) Hx(i,i) += lambda;
                    lambda *= 10;
                    llt.compute(Hx);
                } 

                // Clamp the step size and advance
                DomainType dx = -llt.solve(gx);
                double alpha = 1.0; // maximum step size before hitting a bound
                for (Eigen::Index i = 0; i < gx.size(); i++){
                    if (dx[i] == 0.0) continue;
                    alpha = dx[i] > 0 ? std::min(alpha, (_U[i]-x[i])/dx[i]) : std::min(alpha, (_L[i]-x[i])/dx[i]);
                }
                double fx1 = this->_f(x+alpha*dx);
                while (fx1 > fx + c*alpha*gx.dot(dx)){
                    alpha *= 0.5;
                    fx1 = this->_f(x+alpha*dx);
                }
                x += alpha*dx;
                gx = this->_g(x);
        
                // Convergence check
                if (this->inputConverged(x, dx) || this->gradientConverged(gx) || this->outputConverged(fx, fx1)) return COStatus::Success;
                fx = fx1;
                Hx = _H(x);
            }
            return COStatus::NoConvergence;
        }
    }

    void setFunction(const Function& f, const GFunction& g=nullptr) const noexcept override{
        ConstrainedOptimizer<N>::setFunction(f, g);
        if constexpr (N == 1) _H = [this](const DomainType& x){ return df2(this->_f, x); };
        else _H = [this](const DomainType& x){ return Hessian(this->_f, x); };
    }

    void setFunction(const Function& f, const GFunction& g, const HFunction& H) const noexcept{
        ConstrainedOptimizer<N>::setFunction(f, g);
        _H = H;
    }

    protected:
    mutable DomainType _L; // lower bound
    mutable DomainType _U; // upper bound
    mutable HFunction _H;
};

#endif