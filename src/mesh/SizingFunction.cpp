#include "mesh/SizingFunction.h"

#include <algorithm>
#include <cmath>

namespace mesh {

double SizingFunction::operator()(double d, double s) const
{
    d = std::max(0.0, d);

    // Surface Sizing
    const double gLE = 1 - std::exp(-std::pow((s - sLE_) / wLE_, 2));
    const double gTE = 1 - std::exp(-std::pow((s - sTE_) / wTE_, 2));
    
    // Surface Sizing variation
    const double hSurf = hMin_ + (hMid_ - hMin_) * gLE * gTE;

    // Growth Scaling
    const double gDeltaLE = 1 - std::exp(-std::pow((s - sLE_) / wDelta_, 2));
    const double gDeltaTE = 1 - std::exp(-std::pow((s - sTE_) / wDelta_, 2));

    // Growth Scaling along Surface
    const double delta = deltaMin_ + (deltaMax_ - deltaMin_) * gDeltaLE * gDeltaTE;

    // Final Sizing Function
    const double h = hSurf + (hInf_ - hSurf) * (1 - std::exp(-d / delta));
    return h;
}


} // namespace mesh

