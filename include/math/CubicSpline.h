#ifndef CUBIC_SPLINE_H
#define CUBIC_SPLINE_H

#include "Eigen/Dense"

class CubicSpline{
    public:
    CubicSpline(const Eigen::VectorXd& X, const Eigen::VectorXd& Y, double tol=1e-12);
    
    double sMax() const noexcept{ return _S[_S.size()-1]; }
    Eigen::Vector2d eval(double s) const noexcept;
    Eigen::Vector2d evalDiff(double s) const noexcept;

    //protected:
    Eigen::VectorXd slopes(const Eigen::VectorXd&) noexcept; 
    void splineFit() noexcept;
    Eigen::VectorXd _X, _Y; // dependent variables
    Eigen::VectorXd _S; // independent variable
    Eigen::VectorXd _dX, _dY; // derivative
    double _tol;
};

#endif