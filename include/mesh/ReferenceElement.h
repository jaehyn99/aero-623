#ifndef REFERENCE_ELEMENT_H
#define REFERENCE_ELEMENT_H

#include "Eigen/Dense"
class ReferenceElement{
    public:
    ReferenceElement(int p, int r);

    int p() const noexcept{ return _p; }
    int r() const noexcept{ return _r; }
    int Np() const noexcept{ return (_p+1)*(_p+2)/2; }
    const Eigen::Matrix2Xd& xiL() const noexcept{ return _xiL; }

    const Eigen::Matrix2Xd& intXi() const noexcept{ return _intXi; }
    const Eigen::VectorXd& intW() const noexcept{ return _intW; }
    const Eigen::MatrixXd& intPhi() const noexcept{ return _intPhi; }
    const Eigen::MatrixXd& intPhiXi() const noexcept{ return _intPhiXi; }
    const Eigen::MatrixXd& intPhiEta() const noexcept{ return _intPhiEta; }

    const Eigen::Matrix2Xd& edgeXi(std::size_t i) const noexcept{ return _edgeXi[i]; }
    const Eigen::VectorXd& edgeW() const noexcept{ return _edgeW; }
    const Eigen::MatrixXd& edgePhi(std::size_t i) const noexcept{ return _edgePhi[i]; }
    const Eigen::MatrixXd& edgePhiXi(std::size_t i) const noexcept{ return _edgePhiXi[i]; }
    const Eigen::MatrixXd& edgePhiEta(std::size_t i) const noexcept{ return _edgePhiEta[i]; }

    const Eigen::MatrixXd& M() const noexcept{ return _M; }

    protected:
    int _p; // order of the highest Lagrange basis function used to estimate solution
    int _r; // order of accuracy 
    Eigen::Matrix2Xd _xiL; // Lagrange basis nodes
    
    // Values stored at internal quadrature points
    Eigen::Matrix2Xd _intXi; // internal quadrature points
    Eigen::VectorXd _intW; // internal quadrature weights
    Eigen::MatrixXd _intPhi; // shape functions at internal quadrature points
    Eigen::MatrixXd _intPhiXi; // xi-derivative of shape functions at internal quadrature points
    Eigen::MatrixXd _intPhiEta; // eta-derivative of shape functions at internal quadrature points

    // Values stored at edge quadrature points
    std::array<Eigen::Matrix2Xd, 3> _edgeXi; // edge quadrature points
    Eigen::VectorXd _edgeW; // edge quadrature weights
    std::array<Eigen::MatrixXd, 3> _edgePhi; // shape functions at edge quadrature points
    std::array<Eigen::MatrixXd, 3> _edgePhiXi; // xi-derivative of shape functions at edge quadrature points
    std::array<Eigen::MatrixXd, 3> _edgePhiEta; // eta-derivative of shape functions at edge quadrature points

    Eigen::MatrixXd _M; // mass matrix
};

#endif