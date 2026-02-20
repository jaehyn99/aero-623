#pragma once
#include "numericalFlux.hpp"
#include <Eigen/Dense>
#include <cmath>

class RoeFlux : public numericalFlux {
public:
    Eigen::Vector4d operator()(const Eigen::Vector4d& UL,
                     const Eigen::Vector4d& UR,
                     double gamma, const Eigen::Vector2d& n) const override
    {
	double gm1 = gamma - 1.0;

        double rhoL = UL(0);
        double rhoR = UR(0);

        double uL = UL(1)/rhoL;
        double uR = UR(1)/rhoR;

        double vL = UL(2)/rhoL;
        double vR = UR(2)/rhoR;

	double uNL = uL*n(0) + vL*n(1);
	double uNR = uR*n(0) + vR*n(1);

        double pL = gm1*(UL(3) - 0.5*rhoL*(uL*uL + vL*vL));
        double pR = gm1*(UR(3) - 0.5*rhoR*(uR*uR + vR*vR));

	double HL = (UL(3) + pL)/rhoL;
	double HR = (UR(3) + pR)/rhoR;

	double cL = std::sqrt(gamma*pL/rhoL);
	double cR = std::sqrt(gamma*pR/rhoR);

	Eigen::Matrix<double, 4, 1> FL; Eigen::Matrix<double, 4, 1> FR; 
	Eigen::Matrix<double, 4, 1> F;

	FL(0) = rhoL*uNL;
	FL(1) = UL(1)*uNL + pL*n(0);
	FL(2) = UL(2)*uNL + pL*n(1);
	FL(3) = (UL(3) + pL)*uNL;

	FR(0) = rhoR*uNR;
	FR(1) = UR(1)*uNR + pR*n(0);
	FR(2) = UR(2)*uNR + pR*n(1);
	FR(3) = (UR(3) + pR)*uNR;

	Eigen::Matrix<double, 4, 1> du = UR - UL;

	double di = std::sqrt(rhoR/rhoL);
	double d1 = 1.0/(1.0 + di);

	double ui = (di*uR + uL)*d1;
	double vi = (di*vR + vL)*d1;
	double Hi = (di*HR + HL)*d1;

	double af = 0.5*(ui*ui + vi*vi);
	double ucp = ui*n(0) + vi*n(1);
	double c2 = gm1*(Hi - af);
	double ci = std::sqrt(c2);
	double ci1 = 1.0/ci;

	Eigen::Matrix<double, 3, 1> l;
	l(0) = ucp + ci;
	l(1) = ucp - ci;
	l(2) = ucp;
	
	// double lmax = std::max(l(0), l(1), l(2));
	// double 
	double epsilon = 0.1*ci;

	for (int ii = 0; ii < 3; ii++) {
		if (std::abs(l(ii)) < epsilon) {
        		l(ii) = 0.5 * (epsilon + l(ii)*l(ii)/epsilon);
    		} else {
        		l(ii) = std::abs(l(ii));
    		}
	}

	double s1 = 0.5*(l(0) + l(1));
	double s2 = 0.5*(l(0) - l(1));

	double G1 = gm1*(af*du(0) - ui*du(1) - vi*du(2) + du(3));
	double G2 = -ucp*du(0) + du(1)*n(0) + du(2)*n(1);

	double C1 = G1*(s1 - l(2))*ci1*ci1 + G2*s2*ci1;
	double C2 = G1*s2*ci1 + G2*(s1 - l(2));

	F(0) = 0.5*(FL(0) + FR(0)) - 0.5*(l(2)*du(0) + C1);
	F(1) = 0.5*(FL(1) + FR(1)) - 0.5*(l(2)*du(1) + C1*ui + C2*n(0));
	F(2) = 0.5*(FL(2) + FR(2)) - 0.5*(l(2)*du(2) + C1*vi + C2*n(1));
	F(3) = 0.5*(FL(3) + FR(3)) - 0.5*(l(2)*du(3) + C1*Hi + C2*ucp);


        return F;
    }
};
