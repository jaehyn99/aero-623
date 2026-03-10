#pragma once
#include <Eigen/Dense>

class shapeFunctions {
	public:
		shapeFunctions(int p) : p_(p) {}
		virtual Eigen::MatrixXd getShapeFuncCoeffs(int p_) const;
		virtual Eigen::MatrixXd getShapeFuncXiCoeffs(int p_) const;
		virtual Eigen::MatrixXd getShapeFuncEtaCoeffs(int p_) const;
		virtual Eigen::MatrixXd getRefLagrangePoints(int p_) const;
	protected:
		int p_;
};
