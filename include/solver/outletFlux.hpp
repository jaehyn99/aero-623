#pragma once
#include "boundaryFlux.hpp"
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>

class outletFlux : public boundaryFlux {
private:
    double p0_;

public:
    outletFlux(double p0)
        : p0_(p0) {}

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
        const double pP = std::max(1e-14, gm1 * (rhoEP - 0.5 * rhoP * (uP * uP + vP * vP)));

        const double pBTarget = std::max(1e-14, p0_);
        const double SP = pBTarget / std::pow(rhoP, gamma);
        const double rhoB = std::max(1e-14, std::pow(pBTarget / std::max(1e-14, SP), 1.0 / gamma));
        const double cB = std::sqrt(gamma * pBTarget / rhoB);

        const double uNP = uP * n(0) + vP * n(1);
        const double cP = std::sqrt(gamma * pP / rhoP);
        const double JP = uNP + 2.0 * cP / gm1;
        const double uNB = JP - 2.0 * cB / gm1;

        const double uB = uP - n(0) * (uNP - uNB);
        const double vB = vP - n(1) * (uNP - uNB);
        const double rhoEB = pBTarget / gm1 + 0.5 * rhoB * (uB * uB + vB * vB);

        Eigen::Vector4d F;
        const double MB = std::sqrt(uB * uB + vB * vB) / std::max(1e-14, cB);
        if (MB < 1.0) {
            F(0) = rhoB * (uB * n(0) + vB * n(1));
            F(1) = (rhoB * uB * uB + pBTarget) * n(0) + rhoB * uB * vB * n(1);
            F(2) = rhoB * uB * vB * n(0) + (rhoB * vB * vB + pBTarget) * n(1);
            F(3) = uB * (rhoEB + pBTarget) * n(0) + vB * (rhoEB + pBTarget) * n(1);
        } else {
            F(0) = rhoP * (uP * n(0) + vP * n(1));
            F(1) = (rhoP * uP * uP + pP) * n(0) + rhoP * uP * vP * n(1);
            F(2) = rhoP * uP * vP * n(0) + (rhoP * vP * vP + pP) * n(1);
            F(3) = uP * (rhoEP + pP) * n(0) + vP * (rhoEP + pP) * n(1);
        }

        return F;
    }
};
