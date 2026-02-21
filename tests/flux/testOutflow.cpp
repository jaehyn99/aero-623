#include <iostream>
#include <Eigen/Dense>
#include "solver/outletFlux.hpp"

int main() {
	double pB = 0.1; 
	outletFlux outlet(pB);
	Eigen::Vector4d UP;
	UP << 1.0,   1.0, 0.0, 2.5;
	
	Eigen::Vector2d n;
	n << 1.0, 0.0;

	double gamma = 1.4;

	Eigen::Vector4d FOutlet = outlet(UP, gamma, n);
	
	std::cout << FOutlet << std::endl;
}
