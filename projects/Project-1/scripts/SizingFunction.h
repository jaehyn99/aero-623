#pragma once

#include <algorithm>
#include <cmath>

namespace mesh {

/**
 * h(d, x): target edge length as a function of
 *  - wall distance d >= 0
 *  - blade surface coordinate x (projected x-location on blade)
 *
 * Mathematical form (Python-equivalent):
 *
 *   gLE(x) = 1 - exp(-((x - sLE)/wLE)^2)
 *   gTE(x) = 1 - exp(-((x - sTE)/wTE)^2)
 *
 *   h_surf(x) = hMin + (hMid - hMin) * gLE * gTE
 *
 *   gDLE(x) = 1 - exp(-((x - sLE)/wDelta)^2)
 *   gDTE(x) = 1 - exp(-((x - sTE)/wDelta)^2)
 *
 *   delta(x) = deltaMin + (deltaMax - deltaMin) * gDLE * gDTE
 *
 *   h(d,x) = h_surf(x) + (hInf - h_surf(x)) * (1 - exp(-d / delta(x)))
 *
 * Parameters:
 *   hMin       : minimum wall spacing (>0)
 *   hMid       : surface spacing away from LE/TE (>= hMin)
 *   hInf       : far-field spacing (>= hMid)
 *
 *   sLE, sTE   : leading / trailing edge x-locations
 *
 *   wLE, wTE   : x-widths controlling LE / TE surface refinement
 *   wDelta     : x-width controlling wall-growth scaling
 *
 *   deltaMin   : minimum wall-growth length scale (>0)
 *   deltaMax   : maximum wall-growth length scale (>= deltaMin)
 */
class SizingFunction {
public:
    SizingFunction(double hMin,
                   double hMid,
                   double hInf,
                   double sLE,
                   double sTE,
                   double wLE,
                   double wTE,
                   double wDelta,
                   double deltaMin,
                   double deltaMax)
        : hMin_(hMin), hMid_(hMid), hInf_(hInf),
          sLE_(sLE), sTE_(sTE),
          wLE_(wLE), wTE_(wTE), wDelta_(wDelta),
          deltaMin_(deltaMin), deltaMax_(deltaMax)
    {
        // --- sanitize ---
        if (hMin_ <= 0.0) hMin_ = 1e-12;
        if (hMid_ < hMin_) hMid_ = hMin_;
        if (hInf_ < hMid_) hInf_ = hMid_;

        if (wLE_ <= 0.0)    wLE_    = 1e-12;
        if (wTE_ <= 0.0)    wTE_    = 1e-12;
        if (wDelta_ <= 0.0) wDelta_ = 1e-12;

        if (deltaMin_ <= 0.0) deltaMin_ = 1e-12;
        if (deltaMax_ < deltaMin_) deltaMax_ = deltaMin_;
    }

    /// Full sizing function
    double operator()(double d, double x) const;

    /// Convenience: distance-only (no surface bias)
    double operator()(double d) const { return (*this)(d, 0.0); }

    // --- getters ---
    double hMin() const { return hMin_; }
    double hMid() const { return hMid_; }
    double hInf() const { return hInf_; }

    double sLE() const { return sLE_; }
    double sTE() const { return sTE_; }

    double wLE() const { return wLE_; }
    double wTE() const { return wTE_; }
    double wDelta() const { return wDelta_; }

    double deltaMin() const { return deltaMin_; }
    double deltaMax() const { return deltaMax_; }

private:
    double hMin_, hMid_, hInf_;
    double sLE_, sTE_;
    double wLE_, wTE_, wDelta_;
    double deltaMin_, deltaMax_;
};

} // namespace mesh
