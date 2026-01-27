#include "mesh/SizingFunction.h"

#include <algorithm>
#include <cmath>

namespace mesh {

double SizingFunction::operator()(double d, double xb) const
{
    d = std::max(0.0, d);

    // wall transition
    const double r  = d / std::max(1e-12, d0_);
    const double rp = std::pow(r, std::max(1.0, pWall_));
    const double t  = rp / (1.0 + rp);
    double h = hMin_ + (hMax_ - hMin_) * t;

    // LE/TE bump
    const double sig = std::max(1e-12, xSigma_);
    const double gLE = std::exp(-0.5 * std::pow((xb - xLE_) / sig, 2));
    const double gTE = std::exp(-0.5 * std::pow((xb - xTE_) / sig, 2));
    const double g   = std::max(gLE, gTE);

    const double s = std::clamp(edgeStrength_, 0.0, 1.0);

    // distance gate to prevent far-field LE/TE bias
    const double dEdge = 2.0 * d0_;
    const double q = d / std::max(1e-12, dEdge);
    const double gNear = std::exp(-0.5 * q*q);

    h *= (1.0 - s * g * gNear);

    return std::clamp(h, hMin_, hMax_);
}


} // namespace mesh

