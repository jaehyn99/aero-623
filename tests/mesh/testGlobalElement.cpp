#include "mesh/globalElement.hpp"
#include <iostream>
#include <Eigen/Dense>

int main() {
	int p = 2;
	int q = 4;
	Eigen::MatrixXd V(3, 2);
	V << 1, 1,
	     3, 1,
	     2, 2;

	referenceElement refElem(p, q);
	globalElement globElem(V, refElem);
	std::cout << globElem.J[0] << std::endl;
	std::cout << globElem.detJ[0] << std::endl;
	return 0;
}
