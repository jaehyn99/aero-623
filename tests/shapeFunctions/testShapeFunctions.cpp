#include "ShapeFunctions/ShapeFunctions.hpp"
#include <iostream>
#include <Eigen/Dense>

int main() {
	int p = 3;
	ShapeFunctions shape(p);
	Eigen::MatrixXd phiCoeffs = shape.getShapeFuncCoeffs(p);
	Eigen::MatrixXd phiXiCoeffs = shape.getShapeFuncXiCoeffs(p);
	Eigen::MatrixXd phiEtaCoeffs = shape.getShapeFuncEtaCoeffs(p);
	Eigen::MatrixXd xiL = shape.getRefLagrangePoints(p);
	std::cout << phiCoeffs << std::endl;
	std::cout << phiXiCoeffs << std::endl;
	std::cout << phiEtaCoeffs << std::endl;
	std::cout << xiL << std::endl;
	return 0;
}
