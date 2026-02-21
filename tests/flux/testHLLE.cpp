#include <iostream>
#include <Eigen/Dense>
#include "solver/roe_flux.hpp"
#include "solver/hlle_flux.hpp"

int main() {
	RoeFlux roe;
	HLLEFlux hlle;
	Eigen::Vector4d UL, UR;
	UL << 1.0,   0.0, 0.0, 2.5;
	UR << 0.125, 0.0, 0.0, 2.5;
	
	Eigen::Vector2d n;
	n << 1.0, 0.0;

	double gamma = 1.4;

	Eigen::Vector4d FRoe = roe(UL, UR, gamma, n);
	Eigen::Vector4d FHLLE = hlle(UL, UR, gamma, n);

	std::cout << FRoe << std::endl;
	std::cout << FHLLE << std::endl;
}
