#pragma once
#include "boundaryFlux.hpp"
#include <Eigen/Dense>
#include "FluxResult.hpp"

// DEBUG CHECK =================================
auto check = [](const std::string& name, double val) {
    if (std::isnan(val))
        std::cout << name << " = NaN\n";
    else if (!std::isfinite(val))
        std::cout << name << " = Inf\n";
};
// =============================================

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

    FluxResult operator()(
        const Eigen::Vector4d& UP,
        double gamma,
        const Eigen::Vector2d& n
    ) const override {
	double gm1 = gamma - 1;
        double rhoP = UP(0); // Will make this a function based on the transient flag and time soon
	std::cout << "UP = " << UP.transpose() << "\n";
	
	double uP = UP(1)/rhoP;
	double vP = UP(2)/rhoP;
	double rhoEP = UP(3);
	double pP = gm1*(rhoEP - 0.5*rhoP*(uP*uP + vP*vP));
	double cP = std::sqrt(gamma*pP/rhoP);

	// =============================
	check("rhoP", rhoP);
	check("uP", uP);
	check("vP", vP);
	check("pP", pP);
	check("cP", cP);
	// =============================

	double uNP = uP*n(0) + vP*n(1);
	double pT = rho0_*a0_*a0_/gamma;
	double RTT = pT/rho0_;
	double JP = uNP + 2*cP/gm1;
	double dn = n(0)*std::cos(alpha_) + n(1)*std::sin(alpha_);
	double A = gamma*RTT*dn*dn - 0.5*gm1*JP*JP;
	double B = 4*gamma*RTT*dn/gm1;
	double C = 4*gamma*RTT/(gm1*gm1) - JP*JP;

	// ==================================
	double disc = B*B - 4*A*C;
	if (disc < 0.0)
	{
		std::cout << "\n--- Inlet Quadratic Breakdown ---\n";

		std::cout << "rhoP = " << rhoP << "\n";
		std::cout << "uP   = " << uP << "\n";
		std::cout << "vP   = " << vP << "\n";
		std::cout << "pP   = " << pP << "\n";
		std::cout << "cP   = " << cP << "\n";
		std::cout << "uNP  = " << uNP << "\n";

		std::cout << "pT   = " << pT << "\n";
		std::cout << "RTT  = " << RTT << "\n";
		std::cout << "JP   = " << JP << "\n";
		std::cout << "dn   = " << dn << "\n";

		std::cout << "A    = " << A << "\n";
		std::cout << "B    = " << B << "\n";
		std::cout << "C    = " << C << "\n";
		std::cout << "disc = " << disc << "\n";

		std::cout << "B^2      = " << B*B << "\n";
		std::cout << "4*A*C    = " << 4*A*C << "\n";

		std::cout << "---------------------------------\n";
	}
	check("A", A);
	check("B", B);
	check("C", C);
	check("disc", disc);

	double sqrt_disc = std::sqrt(disc);
	check("sqrt_disc", sqrt_disc);
	// ==================================

	double MB1 = (-B + std::sqrt(B*B - 4*A*C))/(2*A);
	double MB2 = (-B - std::sqrt(B*B - 4*A*C))/(2*A);

	// ==================================
	check("MB1", MB1);
	check("MB2", MB2);
	// ==================================

	double MB;
	if (MB1*MB2 < 0) {
		MB = std::max(MB1, MB2);
	}
	else {
		MB = std::min(std::abs(MB1), std::abs(MB2));
	}

	// ================================
	check("MB", MB);
	// ================================

	double RTB = RTT/(1 + 0.5*gm1*MB*MB);
	double pB = pT*std::pow(RTB/RTT, gamma/gm1);
	double rhoB = pB/RTB;
	double cB = std::sqrt(gamma*pB/rhoB);
	double uNB = JP - 2*cB/gm1;
	double uB = MB*cB*std::cos(alpha_);
	double vB = MB*cB*std::sin(alpha_);
	double rhoEB = pB/gm1 + 0.5*rhoB*(uB*uB + vB*vB);

	// ================================
	check("RTB", RTB);
	check("pB", pB);
	check("rhoB", rhoB);
	check("cB", cB);
	check("uNB", uNB);
	check("uB", uB);
	check("vB", vB);
	check("rhoEB", rhoEB);
	// ================================

	Eigen::Matrix<double, 4, 1> F;
	F(0) = rhoB*(uB*n(0) + vB*n(1));
	F(1) = (rhoB*uB*uB + pB)*n(0) + rhoB*uB*vB*n(1);
	F(2) = rhoB*uB*vB*n(0) + (rhoB*vB*vB + pB)*n(1);
	F(3) = uB*(rhoEB + pB)*n(0) + vB*(rhoEB + pB)*n(1);

	// ---- Wave speed estimate (for CFL only) ----
	double unB = uB*n(0) + vB*n(1);
	double lambdaMax = std::abs(unB) + cB;
	
        return {F, lambdaMax};
    }
};
