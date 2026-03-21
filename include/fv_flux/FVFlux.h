#ifndef FV_FLUX_H
#define FV_FLUX_H

#include <Eigen/Dense>

class FVFlux{
public:
virtual ~FVFlux() = default;
    FVFlux(double gamma): _gamma(gamma) {}
    virtual Eigen::Vector4d computeFlux(const Eigen::Vector4d& UL, const Eigen::Vector4d& UR, const Eigen::Vector2d& n) const = 0;
    virtual double computeWaveSpeed(const Eigen::Vector4d& UL, const Eigen::Vector4d& UR, const Eigen::Vector2d& n, double cL, double cR) const = 0;
    double gamma() const noexcept{ return _gamma; }

    protected:
    double _gamma;
};

#endif