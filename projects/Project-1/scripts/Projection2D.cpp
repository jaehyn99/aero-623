#include "mesh/Projection2D.h"
#include "mesh/BladeGeometry.h"

#include <cmath>
#include <algorithm>

namespace mesh {

static constexpr double gr = 0.6180339887498948482; // (sqrt(5)-1)/2

double Projection2D::objectiveUpper(const BladeGeometry& blade,
                                   double s,
                                   const Eigen::Vector2d& p)
{
    const Eigen::Vector2d r = blade.evalUpper(s);
    return (r - p).squaredNorm();
}

double Projection2D::objectiveLower(const BladeGeometry& blade,
                                   double s,
                                   const Eigen::Vector2d& p)
{
    const Eigen::Vector2d r = blade.evalLower(s);
    return (r - p).squaredNorm();
}


double Projection2D::goldenSectionMinUpper(const BladeGeometry& blade,
                                          double a, double b,
                                          const Eigen::Vector2d& p,
                                          double tol, int maxIter)
{
    double c = b - gr * (b - a);
    double d = a + gr * (b - a);

    double fc = objectiveUpper(blade, c, p);
    double fd = objectiveUpper(blade, d, p);

    for (int it = 0; it < maxIter; ++it) {
        if (std::abs(b - a) < tol) break;

        if (fc < fd) {
            b = d;
            d = c;
            fd = fc;
            c = b - gr * (b - a);
            fc = objectiveUpper(blade, c, p);
        } else {
            a = c;
            c = d;
            fc = fd;
            d = a + gr * (b - a);
            fd = objectiveUpper(blade, d, p);
        }
    }

    return 0.5 * (a + b);
}

double Projection2D::goldenSectionMinLower(const BladeGeometry& blade,
                                          double a, double b,
                                          const Eigen::Vector2d& p,
                                          double tol, int maxIter)
{
    double c = b - gr * (b - a);
    double d = a + gr * (b - a);

    double fc = objectiveLower(blade, c, p);
    double fd = objectiveLower(blade, d, p);

    for (int it = 0; it < maxIter; ++it) {
        if (std::abs(b - a) < tol) break;

        if (fc < fd) {
            b = d;
            d = c;
            fd = fc;
            c = b - gr * (b - a);
            fc = objectiveLower(blade, c, p);
        } else {
            a = c;
            c = d;
            fc = fd;
            d = a + gr * (b - a);
            fd = objectiveLower(blade, d, p);
        }
    }

    return 0.5 * (a + b);
}

ProjectionResult Projection2D::projectToBladeTwo(const BladeGeometry& blade,
                                                const Eigen::Vector2d& p,
                                                double tol,
                                                int maxIter) const
{
    // upper
    ProjectionResult ru;
    ru.curveId = 0;
    ru.s = goldenSectionMinUpper(blade, blade.sUmin(), blade.sUmax(), p, tol, maxIter);
    ru.dist2 = objectiveUpper(blade, ru.s, p);
    ru.dist  = std::sqrt(ru.dist2);
    ru.xProj = blade.evalUpper(ru.s);

    // lower
    ProjectionResult rl;
    rl.curveId = 1;
    rl.s = goldenSectionMinLower(blade, blade.sLmin(), blade.sLmax(), p, tol, maxIter);
    rl.dist2 = objectiveLower(blade, rl.s, p);
    rl.dist  = std::sqrt(rl.dist2);
    rl.xProj = blade.evalLower(rl.s);

    return (ru.dist2 <= rl.dist2) ? ru : rl;
}

ProjectionResult Projection2D::projectToUpper(const BladeGeometry& blade,
                                             const Eigen::Vector2d& p,
                                             double tol, int maxIter) const
{
    ProjectionResult r;
    r.curveId = 0;
    r.s = goldenSectionMinUpper(blade, blade.sUmin(), blade.sUmax(), p, tol, maxIter);
    r.dist2 = objectiveUpper(blade, r.s, p);
    r.dist  = std::sqrt(r.dist2);
    r.xProj = blade.evalUpper(r.s);
    return r;
}

ProjectionResult Projection2D::projectToLower(const BladeGeometry& blade,
                                             const Eigen::Vector2d& p,
                                             double tol, int maxIter) const
{
    ProjectionResult r;
    r.curveId = 1;
    r.s = goldenSectionMinLower(blade, blade.sLmin(), blade.sLmax(), p, tol, maxIter);
    r.dist2 = objectiveLower(blade, r.s, p);
    r.dist  = std::sqrt(r.dist2);
    r.xProj = blade.evalLower(r.s);
    return r;
}

} // namespace mesh
