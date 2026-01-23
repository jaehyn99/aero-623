#include "Constants.h"
#include "CubicSpline.h"
#include "TridiagonalMatrix.h"

#include <algorithm>
#include <iostream>

Eigen::VectorXd slopes(const Eigen::VectorXd& r, const Eigen::VectorXd& _S){
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

    std::cout << "A = " << M.D().transpose() << std::endl;
    std::cout << "B = " << M.L().transpose() << std::endl;
    std::cout << "C = " << M.U().transpose() << std::endl;
    std::cout << "D = " << dr.transpose() << std::endl;

    M.solve(dr);
    return dr;
}

double eval(double s, const Eigen::VectorXd& _X, const Eigen::VectorXd& _S, const Eigen::VectorXd& _dX){
    // Interpolates the spline value at point s
    assert (s >= _S[0] && s <= _S[_S.size()-1]);
    if (s == _S[0]) return _X[0];
    if (s == _S[_S.size()-1]) return _X[_X.size()-1];
    std::size_t ind = std::upper_bound(_S.begin(), _S.end(), s) - _S.begin();
    ind--; // index of the interval where s falls

    double dS = _S[ind+1] - _S[ind];
    double t = (s - _S[ind])/dS;
    double xi0p = _dX[ind]*dS - (_X[ind+1]-_X[ind]);
    double xi1p = _dX[ind+1]*dS - (_X[ind+1]-_X[ind]);
    return (1-t)*_X[ind] + t*_X[ind+1] + t*(1-t)*((1-t)*xi0p - t*xi1p);
}

double evalDiff(double s, const Eigen::VectorXd& _X, const Eigen::VectorXd& _S, const Eigen::VectorXd& _dX){
    // Interpolates the spline derivative at point s
    assert (s >= _S[0] && s <= _S[_S.size()-1]);
    if (s == _S[0]) return _dX[0];
    if (s == _S[_S.size()-1]) return _dX[_dX.size()-1];
    std::size_t ind = std::upper_bound(_S.begin(), _S.end(), s) - _S.begin();
    ind--; // index of the interval where s falls

    double dS = _S[ind+1] - _S[ind];
    double t = (s - _S[ind])/dS;
    double xi0p = _dX[ind]*dS - (_X[ind+1]-_X[ind]);
    double xi1p = _dX[ind+1]*dS - (_X[ind+1]-_X[ind]);
    return (_X[ind+1] - _X[ind] + (1-2*t)*((1-t)*xi0p - t*xi1p) - t*(1-t)*(xi0p+xi1p))/dS;       
}

int main(){
    // Eigen::VectorXd S = Eigen::VectorXd::LinSpaced(6, 0, 2*mconst::pi);
    // Eigen::VectorXd X(6);

    Eigen::VectorXd S{{1, 3, 6, 10}};
    Eigen::VectorXd X{{0, 1, -1, 2}};
    for (int i = 0; i < S.size(); i++) X[i] = std::sin(S[i]);
    auto m = slopes(X, S);
    std::cout << m << std::endl;

    // for (double s: {1, 2, 3, 4, 5, 6, 7, 8, 9, 10}){
    //     double x = eval(s, X, S, m);
    //     double dx = evalDiff(s, X, S, m);
    //     std::cout << "s = " << s << ", x(s) = " << x << ", dx(s) = " << dx << std::endl;
    // }
}