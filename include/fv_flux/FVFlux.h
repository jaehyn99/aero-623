#ifndef FV_FLUX_H
#define FV_FLUX_H

#include <Eigen/Dense>

class FVFlux{
public:
virtual ~FVFlux() = default;
    virtual Eigen::Vector4d computeFlux(const Eigen::Vector4d& UL, const Eigen::Vector4d& UR, const Eigen::Vector2d& n) const = 0;
    virtual double computeWaveSpeed(const Eigen::Vector4d& UL, const Eigen::Vector4d& UR, const Eigen::Vector2d& n, double cL, double cR) const = 0;
};

#endif