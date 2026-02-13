#ifndef CUBIC_SPLINE_H
#define CUBIC_SPLINE_H

#include "Eigen/Dense"

class CubicSpline{
    public:
    CubicSpline(const Eigen::VectorXd& X, const Eigen::VectorXd& Y, double tol=1e-6);
    
    double sMax() const noexcept{ return _S[_S.size()-1]; }
    Eigen::Vector2d eval(double s) const noexcept;
    Eigen::Vector2d evalDeriv(double s) const noexcept;
    Eigen::Vector2d evalDeriv2(double s) const noexcept;
    double projection(const Eigen::Vector2d& p) const;

    //protected:
    Eigen::VectorXd slopes(const Eigen::VectorXd&) noexcept;
    Eigen::VectorXd _X, _Y; // dependent variables
    Eigen::VectorXd _S; // independent variable
    Eigen::VectorXd _dX, _dY; // derivative
    double _tol;
};

#endif