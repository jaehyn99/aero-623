#pragma once

#include <algorithm>
#include <cmath>

namespace mesh {

/**
 * h(d, xb): target edge length as function of wall distance d and projected blade x-location xb.
 *
 * - wall transition: hMin -> hMax as d increases (controlled by d0, pWall)
 * - edge refinement: reduces h near xLE/xTE using Gaussian(s) with width xSigma
 *
 * Parameters:
 *   hMin, hMax  : bounds (hMin>0, hMax>=hMin)
 *   d0          : wall transition length scale (>0)
 *   xLE, xTE    : leading/trailing edge x positions
 *   xSigma      : Gaussian width for edge refinement (>0)
 *   edgeStrength: in [0,1], reduction fraction near edges
 *   pWall       : >=1, sharpness of wall transition
 */
class SizingFunction {
public:
    SizingFunction(double hMin, double hMax,
                   double d0,
                   double xLE, double xTE,
                   double xSigma,
                   double edgeStrength = 0.6,
                   double pWall = 2.0)
        : hMin_(hMin), hMax_(hMax), d0_(d0),
          xLE_(xLE), xTE_(xTE), xSigma_(xSigma),
          edgeStrength_(edgeStrength), pWall_(pWall)
    {
        // sanitize
        if (hMin_ <= 0.0) hMin_ = 1e-12;
        if (hMax_ < hMin_) std::swap(hMax_, hMin_);
        if (d0_ <= 0.0) d0_ = 1e-12;
        if (xSigma_ <= 0.0) xSigma_ = 1e-12;
        edgeStrength_ = std::clamp(edgeStrength_, 0.0, 1.0);
        pWall_ = std::max(1.0, pWall_);
    }

    double operator()(double d, double xb) const;
    double operator()(double d) const { return (*this)(d, 0.0); }

    // getters
    double hMin() const { return hMin_; }
    double hMax() const { return hMax_; }
    double d0()   const { return d0_; }
    double xLE()  const { return xLE_; }
    double xTE()  const { return xTE_; }
    double xSigma() const { return xSigma_; }
    double edgeStrength() const { return edgeStrength_; }
    double pWall() const { return pWall_; }

private:
    double hMin_, hMax_, d0_;
    double xLE_, xTE_, xSigma_;
    double edgeStrength_, pWall_;
};

} // namespace mesh
