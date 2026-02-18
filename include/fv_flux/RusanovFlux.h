#ifndef RUSANOV_FLUX_H
#define RUSANOV_FLUX_H

#include "FVFlux.h"

class RusanovFlux : public FVFlux{
	public:
	RusanovFlux(double gamma): _gamma(gamma) {}
	Eigen::Vector4d computeFlux(const Eigen::Vector4d& UL, const Eigen::Vector4d& UR, const Eigen::Vector2d& n) const override;
	double computeWaveSpeed(const Eigen::Vector4d& UL, const Eigen::Vector4d& UR, const Eigen::Vector2d& n, double cL, double cR) const override;

	protected:
	double _gamma;
};

#endif