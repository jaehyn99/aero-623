#ifndef REFERENCE_ELEMENT_H
#define REFERENCE_ELEMENT_H

#include "Eigen/Dense"
class ReferenceElement{
    public:
    ReferenceElement(int p, int r);

    int p() const noexcept{ return _p; }
    int r() const noexcept{ return _r; }
    int Np() const noexcept{ return (_p+1)*(_p+2)/2; }
    int Nq() const noexcept{ return _wQ.size(); }
    Eigen::Matrix2Xd xiL() const noexcept{ return _xiL; }
    Eigen::Matrix2Xd xiQ() const noexcept{ return _xiQ; }
    Eigen::VectorXd wQ() const noexcept{ return _M; }
    Eigen::MatrixXd phi() const noexcept{ return _phi; }
    Eigen::MatrixXd phiXi() const noexcept{ return _phiXi; }
    Eigen::MatrixXd phiEta() const noexcept{ return _phiEta; }
    Eigen::MatrixXd MRef() const noexcept{ return _M; }

    protected:
    int _p; // order of the highest Lagrange basis function used to estimate solution
    int _r; // order of accuracy 
    Eigen::Matrix2Xd _xiL; // Lagrange basis nodes
    Eigen::Matrix2Xd _xiQ; // quadrature points
    Eigen::VectorXd _wQ; // quadrature weights
    Eigen::MatrixXd _phi; // shape functions at quadrature points
    Eigen::MatrixXd _phiXi; // xi-derivative of shape functions at quadrature points
    Eigen::MatrixXd _phiEta; // eta-derivative of shape functions at quadrature points
    Eigen::MatrixXd _M; // mass matrix
};

#endif