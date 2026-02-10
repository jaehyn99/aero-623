#pragma once
#include "flux.hpp"
#include <cmath>

class RoeFlux : public Flux {
public:
    State operator()(const State& UL,
                     const State& UR,
                     double gamma, const Normal& n) const override
    {
        double rhoL = UL.rho;
        double rhoR = UR.rho;

        double uL = UL.momX / rhoL;
        double uR = UR.momX / rhoR;

        double vL = UL.momY / rhoL;
        double vR = UR.momY / rhoR;

	double uNL = uL*n.X + vL*n.Y;
	double uNR = uR*n.X + vR*n.Y;

        double pL = (gamma - 1.0) * (UL.E - 0.5 * rhoL * (uL * uL + vL * vL));
        double pR = (gamma - 1.0) * (UR.E - 0.5 * rhoR * (uL * uL + vR * vR));

	double rHL = UL.E + pL;
	double rHR = UR.E + pR;

	double HL = rHL / rhoL;
	double HR = rHR / rhoR;

	double cL = std::sqrt(gamma * pL / rhoL);
	double cR = std::sqrt(gamma * pR / rhoR);

        State FL {
            rhoL * uNL,
            UL.momX * uNL + pL*n.X,
            UL.momY * uNL + pL*n.Y,
            (UL.E + pL) * uNL
        };

        State FR {
            rhoR * uNR,
            UR.momX * uNL + pR*n.X,
            UR.momY * uNL + pR*n.Y,
            (UR.E + pR) * uNR
        };

	double du = UR - UL;

	double di = std::sqrt(rR / rL);
	double d1 = 1/(1 + di);

	double ui = (di*uR + uL)*d1;
	double vi = (di*vR + vL)*d1;
	double Hi = (di*HR + HL)*d1;

	double af = 0.5*(ui*ui + vi*vi);
	double ucp = ui*n.X + vi*n.Y;
	double c2 = (gamma - 1)*(Hi - af);
	double ci = std::sqrt(c2);
	double ci1 = 1.0/ci;

	std::vector<double> l = {ucp + ci, ucp - ci, ucp};

        double eta2 = 0.5 * std::sqrt(rhoL * rhoR)
                      / std::pow(std::sqrt(rhoL) + std::sqrt(rhoR), 2);

        double vAvg = (std::sqrt(rhoL) * vL + std::sqrt(rhoR) * vR)
                      / (std::sqrt(rhoL) + std::sqrt(rhoR));

        double EL = UL.E / rhoL;
        double ER = UR.E / rhoR;

        double aL2 = (gamma - 1.0) * ((EL + pL) / rhoL - 0.5 * vL * vL);
        double aR2 = (gamma - 1.0) * ((ER + pR) / rhoR - 0.5 * vR * vR);

        double dAvg = std::sqrt(
            (std::sqrt(rhoL) * aL2 + std::sqrt(rhoR) * aR2)
            / (std::sqrt(rhoL) + std::sqrt(rhoR))
            + eta2 * (vR - vL) * (vR - vL)
        );

        double SL = vAvg - dAvg;
        double SR = vAvg + dAvg;

        if (SL >= 0.0) return FL; // Want to do masking later
        if (SR <= 0.0) return FR;

        return {
            (SR * FL.rho - SL * FR.rho + SL * SR * (UR.rho - UL.rho)) / (SR - SL),
            (SR * FL.mom - SL * FR.mom + SL * SR * (UR.mom - UL.mom)) / (SR - SL),
            (SR * FL.E   - SL * FR.E   + SL * SR * (UR.E   - UL.E))   / (SR - SL)
        };
    }
};
