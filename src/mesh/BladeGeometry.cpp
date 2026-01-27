// BladeGeometry.cpp
// Loads upper and lower point files
// Constructs periodic 2D spline parameterized by arc length s
// Provides eval(s) and deriv(s) - tangent

#include "mesh/BladeGeometry.h"
#include "mesh/TridiagSolver.h"

#include <fstream>
#include <vector>
#include <stdexcept>
#include <cmath>
#include <limits>

namespace mesh {

void BladeGeometry::readXY(const std::string& file, Eigen::VectorXd& x, Eigen::VectorXd& y)
{
    std::ifstream in(file);
    if (!in) throw std::runtime_error("BladeGeometry: cannot open file: " + file);

    std::vector<double> xv, yv;
    double a, b;
    while (in >> a >> b) {
        xv.push_back(a);
        yv.push_back(b);
    }
    if (xv.size() < 2) throw std::runtime_error("BladeGeometry: need >=2 points in: " + file);

    x.resize(static_cast<Eigen::Index>(xv.size()));
    y.resize(static_cast<Eigen::Index>(yv.size()));
    for (Eigen::Index i = 0; i < x.size(); ++i) {
        x(i) = xv[static_cast<size_t>(i)];
        y(i) = yv[static_cast<size_t>(i)];
    }
}

Eigen::VectorXd BladeGeometry::splineFit(const Eigen::VectorXd& X, const Eigen::VectorXd& S)
{
    const Eigen::Index n = X.size();
    if (S.size() != n) throw std::invalid_argument("splineFit: size mismatch");
    if (n < 2) throw std::invalid_argument("splineFit: n<2");

    Eigen::VectorXd A = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd D = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd B = Eigen::VectorXd::Zero(n-1);
    Eigen::VectorXd C = Eigen::VectorXd::Zero(n-1);

    // DS[i] = S[i+1] - S[i]
    Eigen::VectorXd DS = Eigen::VectorXd::Zero(n-1);
    for (Eigen::Index i = 0; i < n-1; ++i) DS(i) = S(i+1) - S(i);

    // Boundary row 0
    A(0) = 2.0 * DS(0);
    C(0) = DS(0);
    D(0) = 3.0 * (X(1) - X(0));

    // Interior rows
    for (Eigen::Index i = 1; i < n-1; ++i) {
        B(i-1) = DS(i-1);
        A(i)   = 2.0 * (DS(i-1) + DS(i));
        C(i)   = DS(i-1);
        D(i)   = 3.0 * ( ( (X(i) - X(i-1)) / DS(i-1) ) * DS(i)
                       + ( (X(i+1) - X(i)) / DS(i)   ) * DS(i-1) );
    }

    // Boundary row n-1
    B(n-2) = DS(n-2);
    A(n-1) = 2.0 * DS(n-2);
    D(n-1) = 3.0 * (X(n-1) - X(n-2));

    return solveTridiag(A, B, C, D);
}

double BladeGeometry::sSimps(double Xi, double Xip1,
                            double Yi, double Yip1,
                            double Si, double Sip1,
                            double dXi, double dXip1,
                            double dYi, double dYip1)
{
    // Direct translation of Python sSimps()
    const double DS = (Sip1 - Si);
    const double slopeX = (Xip1 - Xi) / DS;
    const double slopeY = (Yip1 - Yi) / DS;

    const double xi0P = (dXi   - slopeX) * DS;
    const double xi1P = (dXip1 - slopeX) * DS;
    const double yi0P = (dYi   - slopeY) * DS;
    const double yi1P = (dYip1 - slopeY) * DS;

    const double dxi0 = (Xip1 - Xi) + xi0P;
    const double dyi0 = (Yip1 - Yi) + yi0P;
    const double f0 = std::sqrt(dxi0*dxi0 + dyi0*dyi0);

    const double dxi1 = (Xip1 - Xi) - xi0P/4.0 + xi1P/4.0;
    const double dyi1 = (Yip1 - Yi) - yi0P/4.0 + yi1P/4.0;
    const double f1 = std::sqrt(dxi1*dxi1 + dyi1*dyi1);

    const double dxi2 = (Xip1 - Xi) + xi1P;
    const double dyi2 = (Yip1 - Yi) + yi1P;
    const double f2 = std::sqrt(dxi2*dxi2 + dyi2*dyi2);

    return (f0 + 4.0*f1 + f2) / 6.0;
}

void BladeGeometry::fitCurve(const Eigen::VectorXd& X,
                     const Eigen::VectorXd& Y,
                     double eps,
                     Eigen::VectorXd& S,
                     Eigen::VectorXd& dX,
                     Eigen::VectorXd& dY,
                     CubicSpline1D& xSpline,
                     CubicSpline1D& ySpline)
{
    const Eigen::Index n = X.size();
    S.resize(n);
    S(0) = 0.0;
    for (Eigen::Index i=1; i<n; ++i) {
        const double dx = X(i)-X(i-1);
        const double dy = Y(i)-Y(i-1);
        S(i) = S(i-1) + std::sqrt(dx*dx + dy*dy);
    }

    double splineErr = 1.0;
    Eigen::VectorXd Sold = S;

    while (splineErr > eps) {
        dX = BladeGeometry::splineFit(X, S);
        dY = BladeGeometry::splineFit(Y, S);

        Eigen::VectorXd trueS(n);
        trueS(0) = 0.0;
        splineErr = 0.0;

        for (Eigen::Index i=1; i<n; ++i) {
            trueS(i) = trueS(i-1) + BladeGeometry::sSimps(
                X(i-1), X(i),
                Y(i-1), Y(i),
                S(i-1), S(i),
                dX(i-1), dX(i),
                dY(i-1), dY(i)
            );
            splineErr += std::abs(trueS(i) - S(i));
        }

        S = trueS;
        splineErr /= S(n-1);
    }

    xSpline.set(S, X, dX);
    ySpline.set(S, Y, dY);
}

void BladeGeometry::loadAndFitTwo(const std::string& upperFile,
                                 const std::string& lowerFile,
                                 double eps,
                                double lowerYOffset)
{
    readXY(upperFile, XU_, YU_);
    readXY(lowerFile, XL_, YL_);

    YL_.array() += lowerYOffset;

    fitCurve(XU_, YU_, eps, SU_, dXU_, dYU_, xU_spline_, yU_spline_);
    fitCurve(XL_, YL_, eps, SL_, dXL_, dYL_, xL_spline_, yL_spline_);
}

Eigen::Vector2d BladeGeometry::evalUpper(double s) const {
    return { xU_spline_.eval(s), yU_spline_.eval(s) };
}
Eigen::Vector2d BladeGeometry::evalLower(double s) const {
    return { xL_spline_.eval(s), yL_spline_.eval(s) };
}

} // namespace mesh
