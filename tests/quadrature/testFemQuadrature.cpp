#include "quadrature/femQuadrature.hpp"
#include <iostream>
#include <Eigen/Dense>

int main() {
	int q = 1;
	femQuadrature quad(q);
	Eigen::MatrixXd xi = quad.getQuadXi(q);
	Eigen::VectorXd w = quad.getQuadW(q);
	std::cout << xi << w << std::endl;
	return 0;
}
