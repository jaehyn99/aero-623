#include "shapeFunctions/shapeFunctions.hpp"
#include <iostream>
#include <Eigen/Dense>

Eigen::MatrixXd shapeFunctions::getShapeFuncCoeffs(int p_) const {
	if (p_ == 0) {
		Eigen::MatrixXd phiCoeffs (1, 1);
		phiCoeffs << 1;

		return phiCoeffs;
	}

	if (p_ == 1) {
		Eigen::MatrixXd phiCoeffs (3, 3);
		phiCoeffs << 1, -1, -1,
			     0, 1, 0, 
			     0, 0, 1;

		return phiCoeffs;
	}

	else if (p_ == 2) {
		Eigen::MatrixXd phiCoeffs (6, 6);
		phiCoeffs << 1, -3, -3, 2, 4, 2, 
                             0, 4, 0, -4, -4, 0, 
                             0, -1, 0, 2, 0, 0, 
                             0, 0, 4, 0, -4, -4, 
                             0, 0, 0, 0, 4, 0, 
                             0, 0, -1, 0, 0, 2;

		return phiCoeffs;
	}

	else if (p_ == 3) {
		Eigen::MatrixXd phiCoeffs (10, 10);
		phiCoeffs << 1, -11/2, -11/2, 9, 18, 9, -9/2, -27/2, -27/2, -9/2, 
                             0, 9, 0, -45/2, -45/2, 0, 27/2, 27, 27/2, 0, 
                             0, -9/2, 0, 18, 9/2, 0, -27/2, -27/2, 0, 0, 
                             0, 1, 0, -9/2, 0, 0, 9/2, 0, 0, 0, 
                             0, 0, 9, 0, -45/2, -45/2, 0, 27/2, 27, 27/2, 
                             0, 0, 0, 0, 27, 0, 0, -27, -27, 0, 
                             0, 0, 0, 0, -9/2, 0, 0, 27/2, 0, 0, 
                             0, 0, -9/2, 0, 9/2, 18, 0, 0, -27/2, -27/2,
                             0, 0, 0, 0, -9/2, 0, 0, 0, 27/2, 0,
                             0, 0, 1, 0, 0, -9/2, 0, 0, 0, 9/2;

		return phiCoeffs;
	}

	else {
		throw std::runtime_error("p <= 3");
	}
}

Eigen::MatrixXd shapeFunctions::getShapeFuncXiCoeffs(int p_) const {
	if (p_ == 0) {
		Eigen::MatrixXd phiXiCoeffs (1, 1);
		phiXiCoeffs << 0.0;

		return phiXiCoeffs;
	}
	
	else if (p_ == 1) {
		Eigen::MatrixXd phiXiCoeffs (1, 3);
		phiXiCoeffs << -1.0, 1.0, 0.0;

		return phiXiCoeffs;
	}

	else if (p_ == 2) {
		Eigen::MatrixXd phiXiCoeffs (3, 6);
		phiXiCoeffs << -3.0, 4.0, 4.0, 
			       4, -8, -4, 
			       -1, 4, 0, 
			       0, 0, -4, 
			       0, 0, 4, 
			       0, 0, 0;

		return phiXiCoeffs;
	}

	else if (p_ == 3) {
		Eigen::MatrixXd phiXiCoeffs (6, 10);
		phiXiCoeffs << -11/2, 18, 18, -27/2, -27, -27/2, 
			       9, -45, -45/2, 81/2, 54, 27/2, 
                               -9/2, 36, 9/2, -81/2, -27, 0, 
			       1, -9, 0, 27/2, 0, 0,
                               0, 0, -45/2, 0, 27, 27, 
			       0, 0, 27, 0, -54, -27,
                               0, 0, -9/2, 0, 27, 0, 
			       0, 0, 9/2, 0, 0, -27/2,
                               0, 0, -9/2, 0, 0, 27/2, 
			       0, 0, 0, 0, 0, 0;

	       return phiXiCoeffs;	
	}
	
	else {
		throw std::runtime_error("p <= 3");
	}
}

Eigen::MatrixXd shapeFunctions::getShapeFuncEtaCoeffs(int p_) const {
	if (p_ == 0) {
		Eigen::MatrixXd phiEtaCoeffs (1, 1);
		phiEtaCoeffs << 0.0;

		return phiEtaCoeffs;
	}
	
	else if (p_ == 1) {
		Eigen::MatrixXd phiEtaCoeffs (1, 3);
		phiEtaCoeffs << -1, 0, 1;

		return phiEtaCoeffs;
	}

	else if (p_ == 2) {
		Eigen::MatrixXd phiEtaCoeffs (3, 6);
		phiEtaCoeffs << -3, 4, 4,
			        0, -4, 0, 
				0, 0, 0, 
				4, -4, -8, 
				0, 4, 0, 
				-1, 0, 4;

		return phiEtaCoeffs;
	}

	else if (p_ == 3) {
		Eigen::MatrixXd phiEtaCoeffs (6, 10);
		phiEtaCoeffs << -11/2, 18, 18, -27/2, -27, -27/2, 
			        0, -45/2, 0, 27, 27, 0, 
                                0, 9/2, 0, -27/2, 0, 0, 
				0, 0, 0, 0, 0, 0, 
                                9, -45/2, -45, 27/2, 54, 81/2, 
				0, 27, 0, -27, -54, 0,
                                0, -9/2, 0, 27/2, 0, 0, 
				-9/2, 9/2, 36, 0, -27, -81/2, 
                                0, -9/2, 0, 0, 27, 0, 
				1, 0, -9, 0, 0, 27/2;

	       return phiEtaCoeffs;	
	}
	
	else {
		throw std::runtime_error("p <= 3");
	}
}

Eigen::MatrixXd shapeFunctions::getRefLagrangePoints(int p_) const {
	Eigen::MatrixXd xiL((p_ + 1)*(p_ + 2)/2, 2);
	Eigen::MatrixXd V(3, 2);
	V << 0, 0,
	     1, 0, 
	     0, 1;

	for (int ii = 0; ii < p_ + 1; ii++) {
		for (int jj = 0; jj < p_ + 1 - ii; jj++) {
			for (int kk = 0; kk < 2; kk++){
				xiL(ii*(p_ + 1) - ii*(ii - 1)/2 + jj, kk) = V(0, kk) + ii*(V(2, kk) - V(0, kk))/p_ + jj*(V(1, kk) - V(0, kk))/p_;
			}
		}
	}
	return xiL;
}
