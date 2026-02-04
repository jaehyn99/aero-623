#pragma once

#include <Eigen/Dense>
#include <string>
#include <algorithm>          // std::min/std::max
#include "mesh/CubicSpline1D.h"

namespace mesh {

class BladeGeometry {
public:
    // Two-curve main loader (must be implemented in .cpp)
    void loadAndFitTwo(const std::string& upperFile,
                    const std::string& lowerFile,
                    double eps = 1e-12,
                    double lowerYOffset = 0.0);

    // Backward-compatible wrapper
    void loadAndFit(const std::string& upperFile,
                    const std::string& lowerFile,
                    double eps = 1e-12)
    {
        loadAndFitTwo(upperFile, lowerFile, eps);
    }

    // NOTE: sMin/sMax exist only for convenience; DO NOT use with eval(s).
    double sMin() const { return std::min(sUmin(), sLmin()); }
    double sMax() const { return std::max(sUmax(), sLmax()); }

    // Upper curve API
    double sUmin() const { return SU_.size() ? SU_(0) : 0.0; }
    double sUmax() const { return SU_.size() ? SU_(SU_.size()-1) : 0.0; }
    Eigen::Vector2d evalUpper(double s) const;

    // Lower curve API
    double sLmin() const { return SL_.size() ? SL_(0) : 0.0; }
    double sLmax() const { return SL_.size() ? SL_(SL_.size()-1) : 0.0; }
    Eigen::Vector2d evalLower(double s) const;

private:
    static void readXY(const std::string& file,
                       Eigen::VectorXd& x,
                       Eigen::VectorXd& y);

    // shared helpers
    static Eigen::VectorXd splineFit(const Eigen::VectorXd& X,
                                     const Eigen::VectorXd& S);

    static double sSimps(double Xi, double Xip1,
                         double Yi, double Yip1,
                         double Si, double Sip1,
                         double dXi, double dXip1,
                         double dYi, double dYip1);

    static void fitCurve(const Eigen::VectorXd& X,
                         const Eigen::VectorXd& Y,
                         double eps,
                         Eigen::VectorXd& S,
                         Eigen::VectorXd& dX,
                         Eigen::VectorXd& dY,
                         CubicSpline1D& xSpline,
                         CubicSpline1D& ySpline);

    // per-curve storage
    Eigen::VectorXd XU_, YU_, SU_, dXU_, dYU_;
    Eigen::VectorXd XL_, YL_, SL_, dXL_, dYL_;

    CubicSpline1D xU_spline_, yU_spline_;
    CubicSpline1D xL_spline_, yL_spline_;
};

} // namespace mesh
