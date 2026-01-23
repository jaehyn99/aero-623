#ifndef CONSTRAINED_OPTIMIZER_H
#define CONSTRAINED_OPTIMIZER_H

#include "Eigen/Dense"
#include "Derivative.h"
#include <functional>
#include <type_traits>

enum class COStatus{ Success, InvalidArgument, NoConvergence };

template<int N, typename = std::enable_if_t<(N == Eigen::Dynamic || N >= 1)>>
class ConstrainedOptimizer{
    public:
    using DomainType = std::conditional_t<N == 1, double, Eigen::Matrix<double, N, 1>>;
    using Function = std::function<double(const DomainType&)>;
    using GFunction = std::function<DomainType(const DomainType&)>;

    explicit ConstrainedOptimizer(const Function& func, const GFunction& gfunc=nullptr,
                                  double ftol=1.0e-6, double gtol=1.0e-6, double xtol=1.0e-6, std::size_t maxIter=100):
        _f(func),
        _g(gfunc),
        _ftol(ftol),
        _gtol(gtol),
        _xtol(xtol),
        _maxIter(maxIter)  
    {
        if (!gfunc){
            if constexpr(N == 1) _g = [this](const DomainType& x){ return df(_f, x); };
            else _g = [this](const DomainType& x){ return grad(_f, x); };
        }
    }
    virtual ~ConstrainedOptimizer() = default;

    virtual COStatus minimize(DomainType&) const noexcept = 0;
    virtual void setFunction(const Function& f, const GFunction& g=nullptr) const noexcept{
        _f = f;
        if (g) _g = g;
        else{
            if constexpr(N == 1) _g = [this](const DomainType& x){ return df(_f, x); };
            else _g = [this](const DomainType& x){ return grad(_f, x); };
        }
    }
    void setFTol(double ftol) noexcept{ _ftol = ftol; }
    void setGTol(double gtol) noexcept{ _gtol = gtol; }
    void setXTol(double xtol) noexcept{ _xtol = xtol; }
    void setMaxIter(std::size_t maxIter) noexcept{ _maxIter = maxIter; }

    protected:
    mutable Function _f;
    mutable GFunction _g;
    double _ftol; // Tolerance in output (objective function) step size
    double _gtol; // Tolerance in gradient
    double _xtol; // Tolerance in input step size
    std::size_t _maxIter;

    bool outputConverged(double fx0, double fx1) const noexcept{
        if (fx0 == 0.0) return std::abs(fx1) <= _ftol;
        return std::abs(1.0 - fx1/fx0) <= _ftol;
    }

    bool gradientConverged(const DomainType& gx) const noexcept{
        if constexpr(N == 1) return std::abs(gx) <= _gtol;
        else return gx.squaredNorm() <= _gtol*_gtol;
    }

    bool inputConverged(const DomainType& x, const DomainType& dx) const noexcept{
        if constexpr(N == 1){
            if (x == 0) return std::abs(dx) <= _xtol;
            return std::abs(dx/x) <= _xtol;
        }
        else{
            if (x.squaredNorm() == 0.0) return dx.squaredNorm() <= _xtol*_xtol;
            Eigen::Matrix<double, N, 1> delta(dx);
            for (Eigen::Index i = 0; i < x.size(); i++) delta[i] /= (x[i] == 0.0 ? 1 : x[i]);
            return delta.squaredNorm() <= _xtol*_xtol;
        }
    }
};
#endif