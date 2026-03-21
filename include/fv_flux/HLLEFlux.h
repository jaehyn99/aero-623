#ifndef HLLE_FLUX_H
#define HLLE_FLUX_H

#include "FVFlux.h"

class HLLEFlux: public FVFlux{
        public:
        HLLEFlux(double gamma): FVFlux(gamma) {}
        Eigen::Vector4d computeFlux(const Eigen::Vector4d& UL, const Eigen::Vector4d& UR, const Eigen::Vector2d& n) const override;
        double computeWaveSpeed(const Eigen::Vector4d& UL, const Eigen::Vector4d& UR, const Eigen::Vector2d& n, double cL, double cR) const override;
};
#endif