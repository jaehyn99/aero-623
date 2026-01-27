#pragma once

#include <Eigen/Dense>

namespace mesh {

class BladeGeometry;

struct ProjectionResult {
    int curveId = -1;                    // 0=upper, 1=lower
    double s = 0.0;
    double dist2 = 0.0;
    double dist = 0.0;
    Eigen::Vector2d xProj = Eigen::Vector2d::Zero();
};

class Projection2D {
public:
    ProjectionResult projectToBladeTwo(const BladeGeometry& blade,
                                       const Eigen::Vector2d& p,
                                       double tol = 1e-14,
                                       int maxIter = 200) const;

    // Backward-compatible: returns min-distance among upper/lower
    ProjectionResult projectToBlade(const BladeGeometry& blade,
                                    const Eigen::Vector2d& p,
                                    double tol = 1e-14,
                                    int maxIter = 200) const
    {
        return projectToBladeTwo(blade, p, tol, maxIter);
    }

    ProjectionResult projectToUpper(const BladeGeometry& blade,
                                const Eigen::Vector2d& p,
                                double tol = 1e-14,
                                int maxIter = 200) const;

    ProjectionResult projectToLower(const BladeGeometry& blade,
                                    const Eigen::Vector2d& p,
                                    double tol = 1e-14,
                                    int maxIter = 200) const;

private:
    static double objectiveUpper(const BladeGeometry& blade,
                                 double s,
                                 const Eigen::Vector2d& p);

    static double objectiveLower(const BladeGeometry& blade,
                                 double s,
                                 const Eigen::Vector2d& p);

    static double goldenSectionMinUpper(const BladeGeometry& blade,
                                        double a, double b,
                                        const Eigen::Vector2d& p,
                                        double tol, int maxIter);

    static double goldenSectionMinLower(const BladeGeometry& blade,
                                        double a, double b,
                                        const Eigen::Vector2d& p,
                                        double tol, int maxIter);
                                
};

} // namespace mesh
