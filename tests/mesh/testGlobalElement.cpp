#include "mesh/globalElement.hpp"
#include <iostream>
#include <Eigen/Dense>

int main() {
	int p = 1;
	int q = 2;
	Eigen::MatrixXd V(3, 2);
	V << 1, 1,
	     3, 1,
	     2, 2;

	referenceElement refElem(p, q);
	globalElement globElem(V, refElem);
	std::cout << globElem.M << std::endl;
	return 0;
}
