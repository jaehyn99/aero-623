#pragma once
#include "boundaryFlux.hpp"
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>

class wallFlux : public boundaryFlux {
public:
    Eigen::Vector4d operator()(
        const Eigen::Vector4d& UP,
        double gamma,
        const Eigen::Vector2d& nRaw
    ) const override {
        const double gm1 = gamma - 1.0;
        const double nmag = std::max(1e-14, nRaw.norm());
        const Eigen::Vector2d n = nRaw / nmag;

        const double rhoP = std::max(1e-14, UP(0));
        const double uP = UP(1) / rhoP;
        const double vP = UP(2) / rhoP;
        const double rhoEP = UP(3);

        const double uB = uP - n(0) * (uP * n(0) + vP * n(1));
        const double vB = vP - n(1) * (uP * n(0) + vP * n(1));
        const double pB = std::max(1e-14, gm1 * (rhoEP - 0.5 * rhoP * (uB * uB + vB * vB)));

        Eigen::Vector4d F;
        F(0) = 0.0;
        F(1) = pB * n(0);
        F(2) = pB * n(1);
        F(3) = 0.0;
        return F;
    }
};
