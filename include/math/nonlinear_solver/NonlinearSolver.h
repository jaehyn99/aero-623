#ifndef NONLINEAR_SOLVER_H
#define NONLINEAR_SOLVER_H

#include <Eigen/Dense>
#include <functional>
#include <type_traits>

enum class NLStatus{ Success, InvalidArgument, SingularityError, NoConvergence };

template<int N, int M, typename... Args>
class NonlinearSolver{
    static_assert(N == Eigen::Dynamic || N >= 1, "N must be Eigen::Dynamic or >= 1");
    static_assert(M == Eigen::Dynamic || M >= N, "M must be Eigen::Dynamic or >= N");

    public:
    using DomainType = std::conditional_t<N == 1, double, Eigen::Matrix<double, N, 1>>;
    using RangeType = std::conditional_t<M == 1, double, Eigen::Matrix<double, M, 1>>;
    using Function = std::function<RangeType(const DomainType&)>;

    explicit NonlinearSolver(const Function& func, double ftol=1.0e-6, double xtol=1.0e-6, std::size_t maxIter=100):
        _f(func),
        _ftol(ftol),
        _xtol(xtol),
        _maxIter(maxIter)    
    {}
    virtual ~NonlinearSolver() = default;

    virtual NLStatus solve(DomainType&, Args&&...) const noexcept = 0;
    virtual void setFunction(const Function& f) const noexcept{ _f = f; }
    void setFTol(double ftol) noexcept{ _ftol = ftol; }
    void setXTol(double xtol) noexcept{ _xtol = xtol; }
    void setMaxIter(std::size_t maxIter) noexcept{ _maxIter = maxIter; }

    protected:
    mutable Function _f;
    double _ftol; // Tolerance in output
    double _xtol; // Tolerance in input step size
    std::size_t _maxIter;

    bool outputConverged(const RangeType& fx) const noexcept{
        if constexpr(M == 1) return std::abs(fx) <= _ftol;
        else return fx.squaredNorm() <= _ftol*_ftol;
    }

    bool inputConverged(const DomainType& x, const DomainType& dx) const noexcept{
        if constexpr(N == 1){
            if (x == 0.0) return std::abs(dx) <= _xtol;
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
