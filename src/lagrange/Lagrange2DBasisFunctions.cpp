#include "Lagrange2DBasisFunctions.h"
#include <iostream>

Lagrange2DBasisFunctions::Lagrange2DBasisFunctions(int p):
    _p(p),
    _Np((p+1)*(p+2)/2),
    _phi(Eigen::MatrixXd::Zero(_Np, _Np)),
    _phix(Eigen::MatrixXd::Zero(_Np, p*(p+1)/2)),
    _phiy(Eigen::MatrixXd::Zero(_Np, p*(p+1)/2)), // Empty matrices if p = 0
    _nodes(getLagrangeNodes())
{
    assert(p >= 0 && p <= 3);
	if (_p == 0) _phi << 1;
    else if (_p == 1) {
		// Variable order: 1, x, y
		_phi << 1, -1, -1,
			    0,  1,  0, 
			    0,  0,  1;
	}

	else if (_p == 2) {
		// Variable order: 1, x, y, x^2, xy, y^2
		_phi << 1, -3, -3,  2,  4,  2, 
                0,  4,  0, -4, -4,  0, 
                0, -1,  0,  2 , 0,  0, 
                0,  0,  4,  0, -4, -4, 
                0,  0,  0,  0,  4,  0, 
                0,  0, -1,  0,  0,  2;
	}

	else if (_p == 3) {
		// Variable order: 1, x, y, x^2, xy, y^2, x^3, x^2y, xy^2, y^3
		_phi << 1,  -5.5,   -5.5, 	  9, 	18, 	9,  -4.5, -13.5, -13.5,  -4.5, 
                0, 	   9, 	   0, -22.5, -22.5, 	0,  13.5, 	 27,  13.5, 	0, 
                0,  -4.5, 	   0, 	 18,   4.5, 	0, -13.5, -13.5, 	 0,    	0, 
                0, 	   1, 	   0,  -4.5, 	 0, 	0, 	 4.5, 	  0, 	 0,    	0, 
                0, 	   0, 	   9, 	  0, -22.5, -22.5, 	   0,  13.5, 	27,  13.5, 
                0, 	   0, 	   0, 	  0, 	27, 	0, 	   0, 	-27,   -27,    	0, 
                0, 	   0, 	   0, 	  0,  -4.5, 	0, 	   0,  13.5, 	 0,    	0, 
                0, 	   0, 	-4.5, 	  0,   4.5,    18, 	   0, 	  0, -13.5, -13.5,
                0, 	   0, 	   0, 	  0,  -4.5, 	0, 	   0, 	  0,  13.5, 	0,
                0, 	   0, 	   1, 	  0, 	 0,  -4.5, 	   0, 	  0, 	 0,   4.5;
	}

    if (_p >= 1){
        // Variable order: 1
        _phix.col(0) = _phi.col(1);
        // Variable order: 1
        _phiy.col(0) = _phi.col(2);
    }
    if (_p >= 2){
        // Variable order: (1), 2x, y, first column already updated
        _phix.col(1) = 2*_phi.col(3);
        _phix.col(2) = _phi.col(4);
        // Variable order: (1), x, 2y, first column already updated
        _phiy.col(1) = _phi.col(4);
        _phiy.col(2) = 2*_phi.col(5);
    }
    if (_p == 3){
        // Variable order: (1, 2x, y), 3x^2, 2xy, y^2, first three columns already updated
        _phix.col(3) = 3*_phi.col(6);
        _phix.col(4) = 2*_phi.col(7);
        _phix.col(5) = _phi.col(8);
        // Variable order: (1, x, 2y), x^2, 2xy, 3y^2, first three columns already updated
        _phiy.col(3) = _phi.col(7);
        _phiy.col(4) = 2*_phi.col(8);
        _phiy.col(5) = 3*_phi.col(9);        
    }
}

Eigen::VectorXd Lagrange2DBasisFunctions::evalPhi(double x, double y) const noexcept{
    Eigen::MatrixXd eval = _phi;
    if (_p >= 1){
        eval.col(1) *= x;
        eval.col(2) *= y;
    }
    if (_p >= 2){
        eval.col(3) *= x*x;
        eval.col(4) *= x*y;
        eval.col(5) *= y*y;
    }
    if (_p == 3){
        eval.col(6) *= x*x*x;
        eval.col(7) *= x*x*y;
        eval.col(8) *= x*y*y;
        eval.col(9) *= y*y*y;      
    }
    return eval.rowwise().sum();
}

Eigen::VectorXd Lagrange2DBasisFunctions::evalPhiX(double x, double y) const noexcept{
    Eigen::MatrixXd eval = _phix;
    if (_p >= 2){
        eval.col(1) *= x;
        eval.col(2) *= y;
    }
    if (_p == 3){
        eval.col(3) *= x*x;
        eval.col(4) *= x*y;
        eval.col(5) *= y*y;   
    }
    return eval.rowwise().sum();
}

Eigen::VectorXd Lagrange2DBasisFunctions::evalPhiY(double x, double y) const noexcept{
    Eigen::MatrixXd eval = _phiy;
    if (_p >= 2){
        eval.col(1) *= x;
        eval.col(2) *= y;
    }
    if (_p == 3){
        eval.col(3) *= x*x;
        eval.col(4) *= x*y;
        eval.col(5) *= y*y;   
    }
    return eval.rowwise().sum();
}

Eigen::Matrix2Xd Lagrange2DBasisFunctions::getLagrangeNodes() const noexcept{
    Eigen::Matrix2Xd xiL(2, _Np);
    if (_p == 0) xiL << 0.5, 0.5;
    else{
        Eigen::Matrix<double, 2, 3> V;
        V << 0, 1, 0,
            0, 0, 1;

        for (int ii = 0; ii < _p + 1; ii++) {
            for (int jj = 0; jj < _p + 1 - ii; jj++) {
                for (int kk = 0; kk < 2; kk++){
                    xiL(kk, ii*(_p + 1) - ii*(ii - 1)/2 + jj) = V(kk, 0) + ii*(V(kk, 2) - V(kk, 0))/_p + jj*(V(kk, 1) - V(kk, 0))/_p;
                }
            }
        }
    }
	return xiL;
}

double Lagrange2DBasisFunctions::funcEval(double x, double y, const Eigen::VectorXd& coeff){
    assert(coeff.size() == _Np);
    if (_p == 0) return coeff[0];
    return evalPhi(x,y).dot(coeff);
}

Eigen::VectorXd Lagrange2DBasisFunctions::funcEval(double x, double y, const Eigen::MatrixXd& coeff){
    assert(coeff.cols() == _Np);
    if (_p == 0) return coeff.col(0);
    return coeff*evalPhi(x,y); // result is a Ns-by-1 vector
}

double Lagrange2DBasisFunctions::funcXEval(double x, double y, const Eigen::VectorXd& coeff){
    assert(coeff.size() == _Np);
    if (_p == 0) return 0;
    return evalPhiX(x,y).dot(coeff);
}

Eigen::VectorXd Lagrange2DBasisFunctions::funcXEval(double x, double y, const Eigen::MatrixXd& coeff){
    assert(coeff.cols() == _Np);
    if (_p == 0) return coeff.row(0).transpose();
    return coeff*evalPhiX(x,y); // result is a Ns-by-1 vector
}

double Lagrange2DBasisFunctions::funcYEval(double x, double y, const Eigen::VectorXd& coeff){
    assert(coeff.size() == _Np);
    if (_p == 0) return 0;
    // Reduce coeff to obtain only the coefficients that affect 
    return evalPhiY(x,y).dot(coeff);
}

Eigen::VectorXd Lagrange2DBasisFunctions::funcYEval(double x, double y, const Eigen::MatrixXd& coeff){
    assert(coeff.cols() == _Np);
    if (_p == 0) return coeff.row(0).transpose();
    return coeff*evalPhiY(x,y); // result is a Ns-by-1 vector
}