#pragma once
#include <Eigen/Dense>

class Quadrature1D {
	public:
		Quadrature1D(int q) : q_(q) {}
		virtual Eigen::MatrixXd getQuadXi(int q_) const;
		virtual Eigen::VectorXd getQuadW(int q_) const;
		
	protected:
		int q_;
};
