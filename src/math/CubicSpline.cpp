#include "CubicSpline.h"
#include "GaussLegendre.h"
#include "TridiagonalMatrix.h"

#include <algorithm>
#include <iostream>

CubicSpline::CubicSpline(const Eigen::VectorXd& X, const Eigen::VectorXd& Y, double tol):
    _X(Y),
    _Y(Y),
    _S(X.size()),
    _tol(tol)
{
    assert(X.size() == Y.size());
    //splineFit();
}

Eigen::Vector2d CubicSpline::eval(double s) const noexcept{
    // Interpolates the spline value at point s
    assert(s >= _S[0] && s <= _S[_S.size()-1]);
    std::size_t ind = std::upper_bound(_S.begin(), _S.end(), s) - _S.begin();
    ind--; // index of the interval where s falls

    double dS = _S[ind+1] - _S[ind];
    double t = (s - _S[ind])/dS;
    Eigen::Vector2d dr{(1-t)*_X[ind+1] + t*_X[ind], (1-t)*_Y[ind+1] + t*_Y[ind]};
    Eigen::Vector2d xi0p{_dX[ind]*dS - (_X[ind+1]-_X[ind]), _dY[ind]*dS - (_Y[ind+1]-_Y[ind])};
    Eigen::Vector2d xi1p{_dX[ind+1]*dS - (_X[ind+1]-_X[ind]), _dY[ind+1]*dS - (_Y[ind+1]-_Y[ind])};
    return dr + t*(1-t)*((1-t)*xi0p + t*xi1p);
}

Eigen::Vector2d CubicSpline::evalDiff(double s) const noexcept{
    // Interpolates the spline derivative at point s
    assert(s >= _S[0] && s <= _S[_S.size()-1]);
    std::size_t ind = std::upper_bound(_S.begin(), _S.end(), s) - _S.begin();
    ind--; // index of the interval where s falls

    double dS = _S[ind+1] - _S[ind];
    double t = (s - _S[ind])/dS;
    Eigen::Vector2d dr{_X[ind+1] - _X[ind], _Y[ind+1] - _Y[ind]};
    Eigen::Vector2d xi0p{_dX[ind]*dS - (_X[ind+1]-_X[ind]), _dY[ind]*dS - (_Y[ind+1]-_Y[ind])};
    Eigen::Vector2d xi1p{_dX[ind+1]*dS - (_X[ind+1]-_X[ind]), _dY[ind+1]*dS - (_Y[ind+1]-_Y[ind])};
    return dr + (1-2*t)*((1-t)*xi0p - t*xi1p) - t*(1-t)*(xi0p+xi1p);       
}

Eigen::VectorXd CubicSpline::slopes(const Eigen::VectorXd& r) noexcept{
    // 1D cubic fit on r = r(s)
    Eigen::Index sz = r.size();
    Eigen::VectorXd dS = _S.tail(sz-1) - _S.head(sz-1); // dS has size N-1

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
    M.L(sz-2) = dS[sz-3];
    M.D(sz-1) = 2*dS[sz-3];
    dr[sz-1] = 3*(r[sz-1] - r[sz-2]);

    M.solve(dr);
    return dr;
}

void CubicSpline::splineFit() noexcept{
    Eigen::Index sz = _X.size();
    // Initialize the _S array to be store the arc length on each interval if the spline were linear 
    _S[0] = 0.0;
    for (Eigen::Index i = 1; i < sz; i++){
        double dx = _X[i] - _X[i-1];
        double dy = _Y[i] - _Y[i-1];
        _S[i] = _S[i-1] + std::sqrt(dx*dx + dy*dy);
    }

    double err = 1;
    GaussLegendre<> GK15(7);
    auto arcLength = [this](double s){ return this->evalDiff(s).norm(); };

    std::cout << "S = " << _S.transpose() << std::endl;
    
    while (err > _tol){
        _dX = slopes(_X);
        _dY = slopes(_Y);

        Eigen::VectorXd trueS(sz);
        trueS[0] = 0;
        for (Eigen::Index i = 1; i < sz; i++) trueS[i] = trueS[i-1] + GK15.integrate(arcLength, {_S[i-1], _S[i]});

        err = std::sqrt((trueS - _S).squaredNorm() / trueS.squaredNorm());
        _S = trueS;
        std::cout << "S = " << _S.transpose() << std::endl;
    }
}