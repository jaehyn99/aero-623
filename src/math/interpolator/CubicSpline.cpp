/*#include "CubicSpline.h"
#include "GaussKronrod.h"
#include "ProjectedNewton.h"
#include "TridiagonalMatrix.h"

#include <algorithm>

CubicSpline::CubicSpline(const Eigen::VectorXd& X, const Eigen::VectorXd& Y, double tol):
    _X(X),
    _Y(Y),
    _S(X.size()),
    _tol(tol)
{
    Eigen::Index sz = _X.size();
    assert(Y.size() == sz);

    // Initialize the _S array to be store the arc length on each interval if the spline were linear 
    _S[0] = 0.0;
    for (Eigen::Index i = 1; i < sz; i++){
        double dx = _X[i] - _X[i-1];
        double dy = _Y[i] - _Y[i-1];
        _S[i] = _S[i-1] + std::sqrt(dx*dx + dy*dy);
    }

    double err = 1;
    GaussKronrod<> GK15(7);
    auto arcLength = [this](double s){ return evalDeriv(s).norm(); };
   
    while (err > _tol){
        _dX = slopes(_X);
        _dY = slopes(_Y);

        Eigen::VectorXd trueS(sz);
        trueS[0] = 0;
        for (Eigen::Index i = 1; i < sz; i++) trueS[i] = trueS[i-1] + GK15.integrate(arcLength, {_S[i-1], _S[i]});

        err = std::sqrt((trueS - _S).squaredNorm() / trueS.squaredNorm());
        _S = trueS;
    }
}

Eigen::Vector2d CubicSpline::eval(double s) const noexcept{
    std::size_t sz = _X.size();
    assert (s >= _S[0] && s <= _S[sz-1]);
    if (s == _S[0]) return {_X[0], _Y[0]};
    if (s == _S[sz-1]) return {_X[sz-1], _Y[sz-1]};
    std::size_t ind = std::upper_bound(_S.begin(), _S.end(), s) - _S.begin();
    ind--; // index of the interval where s falls

    double dS = _S[ind+1] - _S[ind];
    double t = (s - _S[ind])/dS;
    Eigen::Vector2d r0{_X[ind], _Y[ind]};
    Eigen::Vector2d r1{_X[ind+1], _Y[ind+1]};
    Eigen::Vector2d dr0{_dX[ind], _dY[ind]};
    dr0 = dr0*dS - (r1-r0);
    Eigen::Vector2d dr1{_dX[ind+1], _dY[ind+1]};
    dr1 = dr1*dS - (r1-r0);
    return (1-t)*r0 + t*r1 + t*(1-t)*((1-t)*dr0 - t*dr1);
}

Eigen::Vector2d CubicSpline::evalDeriv(double s) const noexcept{
    // Interpolates the spline derivative at point s
    std::size_t sz = _X.size();
    assert (s >= _S[0] && s <= _S[sz-1]);
    if (s == _S[0]) return {_dX[0], _dY[0]};
    if (s == _S[sz-1]) return {_dX[sz-1], _dY[sz-1]};
    std::size_t ind = std::upper_bound(_S.begin(), _S.end(), s) - _S.begin();
    ind--; // index of the interval where s falls

    double dS = _S[ind+1] - _S[ind];
    double t = (s - _S[ind])/dS;
    Eigen::Vector2d r0{_X[ind], _Y[ind]};
    Eigen::Vector2d r1{_X[ind+1], _Y[ind+1]};
    Eigen::Vector2d dr0{_dX[ind], _dY[ind]};
    dr0 = dr0*dS - (r1-r0);
    Eigen::Vector2d dr1{_dX[ind+1], _dY[ind+1]};
    dr1 = dr1*dS - (r1-r0);
    return (r1 - r0 + (1-4*t+3*t*t)*dr0 - t*(2-3*t)*dr1)/dS;       
}

Eigen::Vector2d CubicSpline::evalDeriv2(double s) const noexcept{
    // Interpolates the spline derivative at point s
    std::size_t sz = _X.size();
    assert (s >= _S[0] && s <= _S[sz-1]);
    if (s == _S[0] || s == _S[sz-1]) return Eigen::Vector2d::Zero();
    std::size_t ind = std::upper_bound(_S.begin(), _S.end(), s) - _S.begin();
    ind--; // index of the interval where s falls

    double dS = _S[ind+1] - _S[ind];
    double t = (s - _S[ind])/dS;
    Eigen::Vector2d r0{_X[ind], _Y[ind]};
    Eigen::Vector2d r1{_X[ind+1], _Y[ind+1]};
    Eigen::Vector2d dr0{_dX[ind], _dY[ind]};
    dr0 = dr0*dS - (r1-r0);
    Eigen::Vector2d dr1{_dX[ind+1], _dY[ind+1]};
    dr1 = dr1*dS - (r1-r0);
    return ((6*t-4)*dr0 - (6*t-2)*dr1)/(dS*dS);       
}

double CubicSpline::projection(const Eigen::Vector2d& p) const{
    // Find the arclength s such that r(s) is closest to point p
    auto L2   = [this, &p](double s){ return (eval(s)-p).squaredNorm(); };
    auto dL2  = [this, &p](double s){ return 2*(eval(s)-p).dot(evalDeriv(s)); };
    auto d2L2 = [this, &p](double s){ return 2*(evalDeriv(s).squaredNorm() + (eval(s)-p).dot(evalDeriv2(s))); };

    ProjectedNewton<1> newton(L2, 0, _S[_S.size()-1], dL2, d2L2);
    double argmin = _S[_S.size()-1]/2; // initial guess
    COStatus status = newton.minimize(argmin);
    switch (status){
        case COStatus::Success:
            return argmin;
        default:
            throw std::runtime_error("ERROR: No convergence.");
    }
}

Eigen::VectorXd CubicSpline::slopes(const Eigen::VectorXd& r) noexcept{
    // 1D cubic fit on r = r(s)
    Eigen::Index sz = r.size();
    Eigen::VectorXd dS = _S.tail(sz-1) - _S.head(sz-1); // dS has size sz-1

    Eigen::VectorXd dr(sz);
    TridiagonalMatrixXd M(sz);

    // First row, zero curvature boundary
    M.D(0) = 2*dS[0]; // main diagonal, A in the lecture notes 
    M.U(0) = dS[0]; // upper off-diagonal, C in the lecture notes
    dr[0] = 3*(r[1] - r[0]);

    // Middle rows
    M.L().head(sz-2) = dS.tail(sz-2); // lower off-diagonal, B in the lecture notes
    M.D().segment(1, sz-2) = 2*(dS.head(sz-2) + dS.tail(sz-2));
    M.U().tail(sz-2) = dS.head(sz-2);
    for (Eigen::Index i = 1; i < sz-1; i++) dr[i] = 3*((r[i]-r[i-1])/dS[i-1]*dS[i] + (r[i+1]-r[i])/dS[i]*dS[i-1]);

    // Last row, zero curvature boundary
    M.L(sz-2) = dS[sz-2];
    M.D(sz-1) = 2*dS[sz-2];
    dr[sz-1] = 3*(r[sz-1] - r[sz-2]);

    M.solve(dr);
    return dr;
}*/