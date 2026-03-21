#ifndef CUBIC_SPLINE_H
#define CUBIC_SPLINE_H

#include "Eigen/Dense"

class CubicSpline{
    public:
    CubicSpline(const Eigen::VectorXd& X, const Eigen::VectorXd& Y, double tol=1e-6);
    
    Eigen::VectorXd& X() noexcept{ return _X; }
    const Eigen::VectorXd& X() const noexcept{ return _X; }
    Eigen::VectorXd& Y() noexcept{ return _Y; }
    const Eigen::VectorXd& Y() const noexcept{ return _Y; }
    Eigen::VectorXd& S() noexcept{ return _S; }
    const Eigen::VectorXd& S() const noexcept{ return _S; }
    Eigen::VectorXd& dX() noexcept{ return _dX; }
    const Eigen::VectorXd& dX() const noexcept{ return _dX; }
    Eigen::VectorXd& dY() noexcept{ return _dY; }
    const Eigen::VectorXd& dY() const noexcept{ return _dY; }

    double sMax() const noexcept{ return _S[_S.size()-1]; }
    Eigen::Vector2d eval(double s) const noexcept;
    Eigen::Vector2d evalDeriv(double s) const noexcept;
    Eigen::Vector2d evalDeriv2(double s) const noexcept;
    double projection(const Eigen::Vector2d& p) const;

    protected:
    Eigen::VectorXd slopes(const Eigen::VectorXd&) noexcept;
    Eigen::VectorXd _X, _Y; // dependent variables
    Eigen::VectorXd _S; // independent variable
    Eigen::VectorXd _dX, _dY; // derivative
    double _tol;
};

#endif