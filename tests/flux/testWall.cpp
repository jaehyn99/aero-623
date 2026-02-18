#include <iostream>
#include <Eigen/Dense>
#include "solver/wallFlux.hpp"

int main() {
	wallFlux wall;
	Eigen::Vector4d UP;
	UP << 1.0,   1.0, 0.0, 2.5;
	
	Eigen::Vector2d n;
	n << 1.0, 1.0;

	double gamma = 1.4;

	Eigen::Vector4d FWall = wall(UP, gamma, n);
	
	std::cout << FWall << std::endl;
}
