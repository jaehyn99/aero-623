#pragma once
#include "boundaryFlux.hpp"
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>

class inletFlux : public boundaryFlux {
private:
    double rho0_;
    double a0_;
    double alpha_;
    double t_;
    bool   transient_;

public:
    inletFlux(double rho0, double a0, double alpha, double t, bool transient)
        : rho0_(rho0), a0_(a0), alpha_(alpha), t_(t), transient_(transient) {}

    Eigen::Vector4d operator()(
        const Eigen::Vector4d& UP,
        double gamma,
        const Eigen::Vector2d& nRaw
    ) const override {
        (void)t_;
        (void)transient_;

        const double gm1 = gamma - 1.0;
        const double nmag = std::max(1e-14, nRaw.norm());
        const Eigen::Vector2d n = nRaw / nmag;

        const double rhoP = std::max(1e-14, UP(0));
        const double uP = UP(1) / rhoP;
        const double vP = UP(2) / rhoP;
        const double rhoEP = UP(3);
        const double pP = std::max(1e-14, gm1 * (rhoEP - 0.5 * rhoP * (uP * uP + vP * vP)));
        const double cP = std::sqrt(gamma * pP / rhoP);
        const double uNP = uP * n(0) + vP * n(1);

        const double pT = std::max(1e-14, rho0_ * a0_ * a0_ / gamma);
        const double RTT = pT / std::max(1e-14, rho0_);
        const double JP = uNP + 2.0 * cP / gm1;

        const double dn = n(0) * std::cos(alpha_) + n(1) * std::sin(alpha_);
        const double A = gamma * RTT * dn * dn - 0.5 * gm1 * JP * JP;
        const double B = 4.0 * gamma * RTT * dn / gm1;
        const double C = 4.0 * gamma * RTT / (gm1 * gm1) - JP * JP;

        double MB = 0.1;
        const double disc = B * B - 4.0 * A * C;
        if (std::abs(A) > 1e-14 && disc >= 0.0) {
            const double sqrtDisc = std::sqrt(std::max(0.0, disc));
            const double MB1 = (-B + sqrtDisc) / (2.0 * A);
            const double MB2 = (-B - sqrtDisc) / (2.0 * A);
            if (std::isfinite(MB1) && std::isfinite(MB2)) {
                if (MB1 * MB2 < 0.0) {
                    MB = std::max(MB1, MB2);
                } else {
                    MB = std::min(std::abs(MB1), std::abs(MB2));
                }
            }
        }
        if (!std::isfinite(MB)) {
            MB = 0.1;
        }

        const double RTB = RTT / std::max(1e-14, 1.0 + 0.5 * gm1 * MB * MB);
        const double pB = std::max(1e-14, pT * std::pow(std::max(1e-14, RTB / RTT), gamma / gm1));
        const double rhoB = std::max(1e-14, pB / RTB);
        const double cB = std::sqrt(gamma * pB / rhoB);

        const double uB = MB * cB * std::cos(alpha_);
        const double vB = MB * cB * std::sin(alpha_);
        const double rhoEB = pB / gm1 + 0.5 * rhoB * (uB * uB + vB * vB);

        Eigen::Vector4d F;
        F(0) = rhoB * (uB * n(0) + vB * n(1));
        F(1) = (rhoB * uB * uB + pB) * n(0) + rhoB * uB * vB * n(1);
        F(2) = rhoB * uB * vB * n(0) + (rhoB * vB * vB + pB) * n(1);
        F(3) = uB * (rhoEB + pB) * n(0) + vB * (rhoEB + pB) * n(1);

        return F;
    }
};
