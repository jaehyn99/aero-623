#include "RoeFlux.h"
#include <iostream>

Eigen::Vector4d RoeFlux::computeFlux(const Eigen::Vector4d& UL, const Eigen::Vector4d& UR, const Eigen::Vector2d& n) const{
	double gm1 = _gamma - 1.0;

    double rhoL = UL(0);
    double uL = UL(1)/rhoL;
    double vL = UL(2)/rhoL;
    double uNL = uL*n(0) + vL*n(1);
    double pL = gm1*(UL(3) - 0.5*rhoL*(uL*uL + vL*vL));
    double HL = (UL(3) + pL)/rhoL;

    double rhoR = UR(0);
    double uR = UR(1)/rhoR;
    double vR = UR(2)/rhoR;
	double uNR = uR*n(0) + vR*n(1);
    double pR = gm1*(UR(3) - 0.5*rhoR*(uR*uR + vR*vR));
	double HR = (UR(3) + pR)/rhoR;

	Eigen::Vector4d FL{rhoL*uNL, UL(1)*uNL + pL*n(0), UL(2)*uNL + pL*n(1), (UL(3) + pL)*uNL};
    Eigen::Vector4d FR{rhoR*uNR, UR(1)*uNR + pR*n(0), UR(2)*uNR + pR*n(1), (UR(3) + pR)*uNR};

	Eigen::Matrix<double, 4, 1> du = UR - UL;

	double di = std::sqrt(rhoR/rhoL);
	double d1 = 1.0/(1.0 + di); // = std::sqrt(rhoL) / (std::sqrt(rhoR) + std::sqrt(rhoL))
	double ui = (di*uR + uL)*d1;
	double vi = (di*vR + vL)*d1;
	double Hi = (di*HR + HL)*d1;

	double af = 0.5*(ui*ui + vi*vi);
	double ucp = ui*n(0) + vi*n(1);
	double ci = std::sqrt(gm1*(Hi - af));

	Eigen::Vector3d l{ucp + ci, ucp - ci, ucp};
	double epsilon = 0.1*ci;
	for (int ii = 0; ii < 3; ii++) {
		if (std::abs(l(ii)) < epsilon) l(ii) = 0.5 * (epsilon + l(ii)*l(ii)/epsilon);
    	else l(ii) = std::abs(l(ii));
	}

	double s1 = 0.5*(l(0) + l(1));
	double s2 = 0.5*(l(0) - l(1));

	double G1 = gm1*(af*du(0) - ui*du(1) - vi*du(2) + du(3));
	double G2 = -ucp*du(0) + du(1)*n(0) + du(2)*n(1);

	double C1 = G1*(s1 - l(2))/(ci*ci) + G2*s2/ci;
	double C2 = G1*s2/ci + G2*(s1 - l(2));

    Eigen::Vector4d F = 0.5*(FL+FR) - 0.5*l(2)*du;
	F(0) -= 0.5*C1;
	F(1) -= 0.5*(C1*ui + C2*n(0));
	F(2) -= 0.5*(C1*vi + C2*n(1));
	F(3) -= 0.5*(C1*Hi + C2*ucp);
    return F;
}

double RoeFlux::computeWaveSpeed(const Eigen::Vector4d& UL, const Eigen::Vector4d& UR, const Eigen::Vector2d& n, double cL, double cR) const{
	double gm1 = _gamma - 1.0;

    double rhoL = UL(0);
    double uL = UL(1)/rhoL;
    double vL = UL(2)/rhoL;
    double pL = gm1*(UL(3) - 0.5*rhoL*(uL*uL + vL*vL));
    double HL = (UL(3) + pL)/rhoL;

    double rhoR = UR(0);
    double uR = UR(1)/rhoR;
    double vR = UR(2)/rhoR;
    double pR = gm1*(UR(3) - 0.5*rhoR*(uR*uR + vR*vR));
	double HR = (UR(3) + pR)/rhoR;

	double di = std::sqrt(rhoR/rhoL);
	double d1 = 1.0/(1.0 + di); // = std::sqrt(rhoL) / (std::sqrt(rhoR) + std::sqrt(rhoL))
	double ui = (di*uR + uL)*d1;
	double vi = (di*vR + vL)*d1;
	double Hi = (di*HR + HL)*d1;

	double af = 0.5*(ui*ui + vi*vi);
	double ucp = ui*n(0) + vi*n(1);
	double ci = std::sqrt(gm1*(Hi - af));

	Eigen::Vector3d l{ucp + ci, ucp - ci, ucp};
	double epsilon = 0.1*ci;
	for (int ii = 0; ii < 3; ii++) {
		if (std::abs(l(ii)) < epsilon) l(ii) = 0.5 * (epsilon + l(ii)*l(ii)/epsilon);
    	else l(ii) = std::abs(l(ii));
	}
	return l.array().abs().maxCoeff();
}