#include "quadrature/femQuadrature.hpp"
#include <iostream>
#include <Eigen/Dense>

Eigen::MatrixXd femQuadrature::getQuadXi(int q_) const {
	if (q_ == 0 || q_ == 1) {
		Eigen::MatrixXd xi (1, 2);
		xi << 1/3.0, 1/3.0;

		return xi;
	}

	else if (q_ == 2) {
		Eigen::MatrixXd xi (3, 2);
		xi << 0.666666666666667, 0.166666666666667, 
		      0.166666666666667, 0.166666666666667, 
		      0.166666666666667, 0.666666666666667; 

		return xi;
	}

	else if (q_ == 3) {
		Eigen::MatrixXd xi (4, 2);
		xi << 0.333333333333333, 0.333333333333333, 
		      0.600000000000000, 0.200000000000000,
		      0.200000000000000, 0.200000000000000, 
		      0.200000000000000, 0.600000000000000;

		return xi;
	}

	else if (q_ == 4) {
		Eigen::MatrixXd xi (6, 2);
		xi << 0.108103018168070, 0.445948490915965, 
		   0.445948490915965, 0.445948490915965,
		   0.445948490915965, 0.108103018168070, 
		   0.816847572980459, 0.091576213509771,
		   0.091576213509771, 0.091576213509771, 
		   0.091576213509771, 0.816847572980459;

		return xi;
	}

	else {
		throw std::runtime_error("q <= 4");
	}
}

Eigen::VectorXd femQuadrature::getQuadW(int q_) const {
	if (q_ == 0 || q_ == 1) {
		Eigen::VectorXd w(1);
		w << 0.5;

		return w;
	}

	else if (q_ == 2) {
		Eigen::VectorXd w (3);
		w << 0.166666666666666, 0.166666666666666, 0.166666666666666;

		return w;
	}

	else if (q_ == 3) {
		Eigen::VectorXd w (4);
		w << -0.281250000000000, 0.260416666666667, 
		      0.260416666666667, 0.260416666666667; 

		return w;
	}

	else if (q_ == 4) {
		Eigen::VectorXd w (6);
		w <<  0.111690794839005, 0.111690794839005, 
		       0.111690794839005, 0.054975871827661,
		       0.054975871827661, 0.054975871827661;

		return w;
	}

	else {
		throw std::runtime_error("q <= 4");
	}
};
