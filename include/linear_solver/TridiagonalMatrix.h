#ifndef TRIDIAGONAL_MATRIX_H
#define TRIDIAGONAL_MATRIX_H

#include "Eigen/Dense"
#include <type_traits>

template<typename Scalar, int N, 
         typename = std::enable_if_t<(N == Eigen::Dynamic || N > 1)>>
class TridiagonalMatrix{
    public:
    using OffDiagonal = Eigen::Matrix<Scalar, (N == Eigen::Dynamic) ? Eigen::Dynamic : N-1, 1>;
    using Diagonal = Eigen::Matrix<Scalar, N, 1>;
    TridiagonalMatrix() = default;
    TridiagonalMatrix(const OffDiagonal& l, const Diagonal& d, const OffDiagonal& u):
        _l(l), _d(d), _u(u)
    {
        if constexpr (N == Eigen::Dynamic) assert(l.size() == d.size()-1 && u.size() == d.size()-1);
    }

    template<int M = N, typename = std::enable_if_t<M == Eigen::Dynamic>>
    TridiagonalMatrix(std::size_t n): _l(n-1), _d(n), _u(n-1)
    {
        assert(n > 1);
    }

    OffDiagonal& L() noexcept{ return _l; }
    const OffDiagonal& L() const noexcept{ return _l; }
    Scalar& L(Eigen::Index i) noexcept{ return _l[i]; }
    const Scalar& L(Eigen::Index i) const noexcept{ return _l[i]; }

    Diagonal& D() noexcept{ return _d; }
    const Diagonal& D() const noexcept{ return _d; }
    Scalar& D(Eigen::Index i) noexcept{ return _d[i]; }
    const Scalar& D(Eigen::Index i) const noexcept{ return _d[i]; }

    OffDiagonal& U() noexcept{ return _u; }
    const OffDiagonal& U() const noexcept{ return _u; }
    Scalar& U(Eigen::Index i) noexcept{ return _u[i]; }
    const Scalar& U(Eigen::Index i) const noexcept{ return _u[i]; }

    template<typename = std::enable_if_t<N == Eigen::Dynamic>>
    void resize(Eigen::Index n) const noexcept{
        assert(n > 1);
        _l.resize(n-1);
        _d.resize(n);
        _u.resize(n-1);
    }

    void solve(Diagonal& rhs) const noexcept{
        assert(rhs.size() == _d.size());
        Diagonal d(_d); // temp copy of the diagonal
        for (Eigen::Index i = 1; i < d.size(); i++){
            Scalar w = _l[i-1]/d[i-1]; 
            d[i] -= w*_u[i-1];
            rhs[i] -= w*rhs[i-1];
        }

        rhs[rhs.size()-1] /= d[d.size()-1];
        for (Eigen::Index i = d.size()-2; i >= 0; i--){
            rhs[i] -= _u[i]*rhs[i+1];
            rhs[i] /= d[i];
        }
    }

    protected:
    OffDiagonal _l;
    Diagonal _d;
    OffDiagonal _u;
};

using TridiagonalMatrix2d = TridiagonalMatrix<double, 2>;
using TridiagonalMatrix3d = TridiagonalMatrix<double, 3>;
using TridiagonalMatrix4d = TridiagonalMatrix<double, 4>;
using TridiagonalMatrixXd = TridiagonalMatrix<double, Eigen::Dynamic>;

#endif