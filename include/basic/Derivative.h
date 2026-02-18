#ifndef DERIVATIVE_H
#define DERIVATIVE_H

#include "Eigen/Dense"
#include <cmath>
#include <limits>
#include <iostream>

// First derivative
// Derivative of a single-variable function
template<typename Callable, typename Scalar,
         typename = std::enable_if_t<std::is_floating_point<Scalar>::value>>
auto df(Callable&& f, Scalar x) -> decltype(f(x)) {
    const auto eps = std::numeric_limits<Scalar>::epsilon();
    Scalar dx = std::cbrt(3*eps) * std::max(std::abs(x), 1.0);
    return (f(x+dx) - f(x-dx))/(2*dx);
}

// Partial derivative of a multi-variable function, with respect to variable x[ind]
template<typename Callable, typename Scalar, int Rows, int Cols,
         typename = std::enable_if_t<std::is_floating_point<Scalar>::value>>
auto df(Callable&& f, const Eigen::Matrix<Scalar, Rows, Cols>& x, Eigen::Index ind) -> decltype(f(x)){
    static_assert(Rows == 1 || Cols == 1, "Parameter x must be a row or column vector.");
    if (ind >= x.size()) throw std::invalid_argument("ERROR: Index out of bound.");

    const auto eps = std::numeric_limits<Scalar>::epsilon();
    Scalar dx = std::cbrt(3*eps) * std::max(std::abs(x[ind]), 1.0);
    auto xp = x; xp[ind] += dx;
    auto xm = x; xm[ind] -= dx;
    return (f(xp) - f(xm))/(2*dx);
}

// First derivative operators
// Gradient
template<typename Callable, typename Scalar, int Rows, int Cols,
         typename = std::enable_if_t<std::is_floating_point<Scalar>::value>>
auto grad(Callable&& f, const Eigen::Matrix<Scalar, Rows, Cols>& x){
    static_assert(Rows == 1 || Cols == 1, "Parameter x must be a row or column vector.");
    static_assert(std::is_floating_point<decltype(f(x))>::value, "Function must return a scalar.");
    
    std::decay_t<decltype(x)> gradf;
    if (Rows == Eigen::Dynamic || Cols == Eigen::Dynamic) gradf.resize(x.size());
    for (Eigen::Index i = 0; i < x.size(); i++) gradf[i] = df(f, x, i);
    return gradf;
}

// Divergent
template<typename Callable, typename Scalar, int Rows, int Cols,
         typename = std::enable_if_t<std::is_floating_point<Scalar>::value>>
Scalar div(Callable&& f, const Eigen::Matrix<Scalar, Rows, Cols>& x){
    static_assert(Rows == 1 || Cols == 1, "Parameter x must be a row or column vector.");
    if (f(x).size() != x.size()) throw std::invalid_argument("ERROR: The domain and range dimensions must match.");

    Scalar divf = 0.0;
    for (Eigen::Index i = 0; i < x.size(); i++) divf += df(f, x, i)[i];
    return divf;
}

// Curl
template<typename Callable, typename Scalar, int Rows, int Cols,
         typename = std::enable_if_t<std::is_floating_point<Scalar>::value>>
auto curl(Callable&& f, const Eigen::Matrix<Scalar, Rows, Cols>& x){
    static_assert(Rows == 1 || Cols == 1, "Parameter x must be a row or column vector.");
    Eigen::Index domain = x.size();
    Eigen::Index range = f(x).size();
    if (domain != 2 && domain != 3)
        throw std::invalid_argument("ERROR: Curl is only defined for the R^3 domain.");
    if (range != 2 && range != 3)
        throw std::invalid_argument("ERROR: Curl is only defined for the R^3 range.");

    std::decay_t<decltype(x)> curlf;
    if (Rows == Eigen::Dynamic || Cols == Eigen::Dynamic) curlf.resize(3);
    curlf.setZero();

    auto df0 = df(f, x, 0);
    auto df1 = df(f, x, 1);
    if (domain == 3){
        auto df2 = df(f, x, 2);
        curlf[0] -= df2[1];
        curlf[1] += df2[0];
    }
    if (range == 3){
        curlf[0] += df1[2];
        curlf[1] -= df0[2];
    }
    curlf[2] = df0[1] - df1[0];
    return curlf;
}

// Jacobian
template<typename Callable, typename Scalar, int Rows, int Cols,
         typename = std::enable_if_t<std::is_floating_point<Scalar>::value>>
auto Jacobian(Callable&& f, const Eigen::Matrix<Scalar, Rows, Cols>& x){
    static_assert(Rows == 1 || Cols == 1, "Parameter x must be a row or column vector.");
    Eigen::Index domain = x.size();
    Eigen::Index range = f(x).size();

    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Jf;
    Jf.resize(range, domain);
    for (Eigen::Index i = 0; i < domain; i++){
        Jf.col(i) = Cols == 1 ? df(f, x, i) : df(f, x, i).transpose();
    }
    return Jf;
}

// Second derivative
// Derivative of a single-variable function
template<typename Callable, typename Scalar,
         typename = std::enable_if_t<std::is_floating_point<Scalar>::value>>
auto df2(Callable&& f, Scalar x) -> decltype(f(x)){
    const auto eps = std::numeric_limits<Scalar>::epsilon();
    Scalar dx = std::sqrt(std::sqrt(48*eps)) * std::max(std::abs(x), 1.0);
    return (f(x+dx) - 2*f(x) + f(x-dx))/(dx*dx);
}

// Partial derivative of a multi-variable function, with respect to variable x[ind]
template<typename Callable, typename Scalar, int Rows, int Cols,
         typename = std::enable_if_t<std::is_floating_point<Scalar>::value>>
auto df2(Callable&& f, const Eigen::Matrix<Scalar, Rows, Cols>& x, Eigen::Index ind1, Eigen::Index ind2) -> decltype(f(x)){
    static_assert(Rows == 1 || Cols == 1, "Parameter x must be a row or column vector.");
    if (ind1 >= x.size() || ind2 >= x.size()) throw std::invalid_argument("ERROR: Index out of bound.");
    
    const auto eps = std::numeric_limits<Scalar>::epsilon();
    Scalar dx = std::sqrt(std::sqrt(48*eps));
    if (ind1 == ind2){ // Own partial derivatives
        dx *= std::max(std::abs(x[ind1]), 1.0);
        auto xp = x; xp[ind1] += dx;
        auto xm = x; xm[ind1] -= dx;
        return (f(xp) - 2*f(x) + f(xm))/(dx*dx);
    }
    // Mixed partial derivatives
    dx /= std::sqrt(std::sqrt(2));
    Scalar dx1 = dx*std::max(std::abs(x[ind1]), 1.0);
    Scalar dx2 = dx*std::max(std::abs(x[ind2]), 1.0);
    auto xpp = x; xpp[ind1] += dx1; xpp[ind2] += dx2;
    auto xpm = x; xpm[ind1] += dx1; xpm[ind2] -= dx2;
    auto xmp = x; xmp[ind1] -= dx1; xmp[ind2] += dx2;
    auto xmm = x; xmm[ind1] -= dx1; xmm[ind2] -= dx2;
    return (f(xpp) - f(xpm) - f(xmp) + f(xmm))/(4*dx*dx);
}

// Second derivative operators
// Scalar Laplacian
template<typename Callable, typename Scalar, int Rows, int Cols,
         typename = std::enable_if_t<std::is_floating_point<Scalar>::value>>
Scalar Laplacian(Callable&& f, const Eigen::Matrix<Scalar, Rows, Cols>& x){
    static_assert(Rows == 1 || Cols == 1, "Parameter x must be a row or column vector.");
    static_assert(std::is_floating_point<decltype(f(x))>::value, "Function must return a scalar.");

    Scalar laplace = 0.0;
    for (Eigen::Index i = 0; i < x.size(); i++){
        auto df2_xi = df2(f, x, i, i);
        laplace += df2_xi;
    }
    return laplace;
}

// Hessian
template<typename Callable, typename Scalar, int Rows, int Cols,
         typename = std::enable_if_t<std::is_floating_point<Scalar>::value>>
auto Hessian(Callable&& f, const Eigen::Matrix<Scalar, Rows, Cols>& x){
    static_assert(Rows == 1 || Cols == 1, "Parameter x must be a row or column vector.");
    static_assert(std::is_floating_point<decltype(f(x))>::value, "Function must return a scalar.");

    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Hf;
    const Eigen::Index n = x.size();
    Hf.resize(n, n);
    for (Eigen::Index i = 0; i < n; i++){
        Hf(i,i) = df2(f, x, i, i);
        for (Eigen::Index j = i+1; j < n; j++){
            Hf(i,j) = df2(f, x, i, j);
            Hf(j,i) = Hf(i,j);
        }
    }
    return Hf;
}

#endif