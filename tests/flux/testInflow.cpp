#include <iostream>
#include <Eigen/Dense>
#include "solver/roeFlux.hpp"
#include "solver/hlleFlux.hpp"
#include "solver/inletFlux.hpp"
#include "solver/outletFlux.hpp"
#include "solver/wallFlux.hpp"

int main() {
	double rho0 = 1.0; 
	double a0 = 2.0;
	double alpha = 0.25;
	double t = 0;
	bool transient = 0;
	inletFlux inlet(rho0, a0, alpha, t, transient);
	Eigen::Vector4d UP;
	UP << 1.0,   1.0, 0.0, 2.5;
	
	Eigen::Vector2d n;
	n << 1.0, 0.0;

	double gamma = 1.4;

	Eigen::Vector4d FInlet = inlet(UP, gamma, n);
	
	std::cout << FInlet << std::endl;
}
