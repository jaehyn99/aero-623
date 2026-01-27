#include "mesh/CubicSpline1D.h"

#include <algorithm>
#include <stdexcept>
#include <cmath>

namespace mesh {

void CubicSpline1D::set(const Eigen::VectorXd& S,
                        const Eigen::VectorXd& X,
                        const Eigen::VectorXd& dX)
{
    if (S.size() < 2) throw std::invalid_argument("CubicSpline1D::set: need >=2 knots");
    if (X.size() != S.size() || dX.size() != S.size())
        throw std::invalid_argument("CubicSpline1D::set: size mismatch");
    S_ = S;
    X_ = X;
    dX_ = dX;
}

Eigen::Index CubicSpline1D::findInterval(double s) const
{
    // clamp to valid range
    if (s <= S_(0)) return 0;
    if (s >= S_(S_.size()-1)) return S_.size()-2;

    // binary search on knots
    Eigen::Index lo = 0;
    Eigen::Index hi = S_.size(); // one past end
    while (lo < hi) {
        Eigen::Index mid = (lo + hi) / 2;
        if (S_(mid) < s) lo = mid + 1;
        else hi = mid;
    }
    Eigen::Index ind = lo - 1;
    if (ind < 0) ind = 0;
    if (ind > S_.size()-2) ind = S_.size()-2;
    return ind;
}

double CubicSpline1D::eval(double s) const
{
    const Eigen::Index i = findInterval(s);
    const double Si  = S_(i);
    const double Sip = S_(i+1);
    const double DSi = Sip - Si;
    const double t   = (s - Si) / DSi;

    const double Xi  = X_(i);
    const double Xip = X_(i+1);

    // Python:
    // xi0P = (dXi[i] - (Xip - Xi)/DSi)*DSi
    // xi1P = (dXi[i+1] - (Xip - Xi)/DSi)*DSi
    const double slope = (Xip - Xi) / DSi;
    const double xi0P  = (dX_(i)   - slope) * DSi;
    const double xi1P  = (dX_(i+1) - slope) * DSi;

    // x = (1 - t)*Xi + t*Xip + (t - t^2)*((1 - t)*xi0P - t*xi1P)
    const double x = (1.0 - t)*Xi + t*Xip
                   + (t - t*t) * ((1.0 - t)*xi0P - t*xi1P);

    return x;
}

double CubicSpline1D::deriv(double s) const
{
    const Eigen::Index i = findInterval(s);
    const double Si  = S_(i);
    const double Sip = S_(i+1);
    const double DSi = Sip - Si;
    const double t   = (s - Si) / DSi;

    const double Xi  = X_(i);
    const double Xip = X_(i+1);

    const double slope = (Xip - Xi) / DSi;
    const double xi0P  = (dX_(i)   - slope) * DSi;
    const double xi1P  = (dX_(i+1) - slope) * DSi;

    // Python dx (but note python returns dx/dt form; then /DSi for dx/ds)
    // dx = Xip - Xi + (1 - 2t)*((1 - t)*xi0P - t*xi1P) - (t - t^2)*(xi0P + xi1P)
    const double dx_dt = (Xip - Xi)
        + (1.0 - 2.0*t) * ((1.0 - t)*xi0P - t*xi1P)
        - (t - t*t) * (xi0P + xi1P);

    return dx_dt / DSi;
}

} // namespace mesh
