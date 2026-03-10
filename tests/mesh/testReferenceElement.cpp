#include "mesh/referenceElement.hpp"
#include <iostream>
#include <Eigen/Dense>

int main() {
	int p = 1;
	int q = 2;
	referenceElement refElem(p, q);
	std::cout << refElem.phiEta << std::endl;
	return 0;
}
