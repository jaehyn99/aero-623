#include "Lagrange1DBasisFunctions.h"

Lagrange1DBasisFunctions::Lagrange1DBasisFunctions(int p):
    _p(p),
    _phi(Eigen::MatrixXd::Zero(_p+1, _p+1)),
    _phix(Eigen::MatrixXd::Zero(_p+1, _p))
{
    assert(p >= 0 && p <= 3);
	if (_p == 0) _phi << 1;
    else if (_p == 1) {
		// Variable order: 1, x
		_phi << 1, -1,
			    0,  1;
	}

	else if (_p == 2) {
		// Variable order: 1, x, x^2
		_phi << 1, -3,  2, 
                0,  4, -4, 
                0, -1,  2;
	}

	else if (_p == 3) {
		// Variable order: 1, x, x^2, x^3,
		_phi << 1,  -5.5,      9,  -4.5,
                0, 	   9,  -22.5,  13.5,
                0,  -4.5, 	  18, -13.5,
                0, 	   1, 	-4.5,   4.5;
	}

    _phix = _phi.rightCols(_p);
    if (_p >= 2) _phix.col(1) *= 2;
    if (_p == 3) _phix.col(2) *= 3;
}

Eigen::VectorXd Lagrange1DBasisFunctions::evalPhi(double x) const noexcept{
    if (_p == 0) return _phi; // 1

    Eigen::MatrixXd eval = _phi;
    if (_p >= 1) eval.col(1) *= x;
    if (_p >= 2) eval.col(2) *= x*x;
    if (_p == 3) eval.col(3) *= x*x*x;
    return eval.rowwise().sum();
}

Eigen::VectorXd Lagrange1DBasisFunctions::evalPhiX(double x) const noexcept{
    if (_p == 0) return Eigen::MatrixXd(); // empty matrix

    Eigen::MatrixXd eval = _phix;
    if (_p >= 2) eval.col(1) *= x;
    if (_p == 3) eval.col(2) *= x*x;
    return eval.rowwise().sum();
}

double Lagrange1DBasisFunctions::funcEval(double x, const Eigen::VectorXd& coeff){
    assert(coeff.size() == _p+1);
    if (_p == 0) return coeff[0];
    return evalPhi(x).dot(coeff);
}

Eigen::VectorXd Lagrange1DBasisFunctions::funcEval(double x, const Eigen::MatrixXd& coeff){
    assert(coeff.cols() == _p+1);
    if (_p == 0) return coeff.col(0);
    return coeff*evalPhi(x); // result is a Ns-by-1 vector
}

double Lagrange1DBasisFunctions::funcXEval(double x, const Eigen::VectorXd& coeff){
    assert(coeff.size() == _p+1);
    if (_p == 0) return 0;
    return evalPhiX(x).dot(coeff);
}

Eigen::VectorXd Lagrange1DBasisFunctions::funcXEval(double x, const Eigen::MatrixXd& coeff){
    assert(coeff.cols() == _p+1);
    if (_p == 0) return coeff.row(0).transpose();
    return coeff*evalPhiX(x); // result is a Ns-by-1 vector
}