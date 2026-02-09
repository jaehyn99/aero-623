#pragma once
#include "flux.hpp"
#include <cmath>

class HLLEFlux : public Flux {
public:
    State operator()(const State& uL,
                     const State& uR,
                     double gamma) const override
    {
        double rhoL = uL.rho;
        double rhoR = uR.rho;

        double vL = uL.mom / rhoL;
        double vR = uR.mom / rhoR;

        double pL = (gamma - 1.0) * (uL.E - 0.5 * rhoL * vL * vL);
        double pR = (gamma - 1.0) * (uR.E - 0.5 * rhoR * vR * vR);

        State FL {
            rhoL * vL,
            rhoL * vL * vL + pL,
            (uL.E + pL) * vL
        };

        State FR {
            rhoR * vR,
            rhoR * vR * vR + pR,
            (uR.E + pR) * vR
        };

        double eta2 = 0.5 * std::sqrt(rhoL * rhoR)
                      / std::pow(std::sqrt(rhoL) + std::sqrt(rhoR), 2);

        double vAvg = (std::sqrt(rhoL) * vL + std::sqrt(rhoR) * vR)
                      / (std::sqrt(rhoL) + std::sqrt(rhoR));

        double EL = uL.E / rhoL;
        double ER = uR.E / rhoR;

        double aL2 = (gamma - 1.0) * ((EL + pL) / rhoL - 0.5 * vL * vL);
        double aR2 = (gamma - 1.0) * ((ER + pR) / rhoR - 0.5 * vR * vR);

        double dAvg = std::sqrt(
            (std::sqrt(rhoL) * aL2 + std::sqrt(rhoR) * aR2)
            / (std::sqrt(rhoL) + std::sqrt(rhoR))
            + eta2 * (vR - vL) * (vR - vL)
        );

        double SL = vAvg - dAvg;
        double SR = vAvg + dAvg;

        if (SL >= 0.0) return FL;
        if (SR <= 0.0) return FR;

        return {
            (SR * FL.rho - SL * FR.rho + SL * SR * (uR.rho - uL.rho)) / (SR - SL),
            (SR * FL.mom - SL * FR.mom + SL * SR * (uR.mom - uL.mom)) / (SR - SL),
            (SR * FL.E   - SL * FR.E   + SL * SR * (uR.E   - uL.E))   / (SR - SL)
        };
    }
};
