// The set of 2D Lagrange basis functions defined on a unit right triangle

#ifndef LAGRANGE_BASIS_FUNCTIONS_H
#define LAGRANGE_BASIS_FUNCTIONS_H

#include "Eigen/Dense"
class LagrangeBasisFunctions{
    public:
    LagrangeBasisFunctions(int p);

    int p() const noexcept{ return _p; }
    int Np() const noexcept{ return _Np; } 
    Eigen::MatrixXd phi() const noexcept{ return _phi; }
    Eigen::MatrixXd phix() const noexcept{ return _phix; }
    Eigen::MatrixXd phiy() const noexcept{ return _phiy; }
    Eigen::Matrix2Xd nodes() const noexcept{ return _nodes; }

    // Evaluation of polynomial and its derivatives
    double funcEval (double x, double y, const Eigen::VectorXd& coeff);
    double funcXEval(double x, double y, const Eigen::VectorXd& coeff);
    double funcYEval(double x, double y, const Eigen::VectorXd& coeff);

    // Multi-D overloads that returns an Ns-by-1 vector
    // coeff is of dimension Ns*Np, where Ns = number of states, Np = number of basis functions;
    Eigen::VectorXd funcEval (double x, double y, const Eigen::MatrixXd& coeff);
    Eigen::VectorXd funcXEval(double x, double y, const Eigen::MatrixXd& coeff);
    Eigen::VectorXd funcYEval(double x, double y, const Eigen::MatrixXd& coeff);

    protected:
    int _p; // polynomial order
    int _Np; // number of basis polynomials
    Eigen::MatrixXd _phi; // matrix of basis function coefficients, each row is a basis with the monomial coefficients given in the columns
    Eigen::MatrixXd _phix; // matrix of x-derivative coefficients, each row is a basis with the monomial coefficients given in the columns
    Eigen::MatrixXd _phiy; // matrix of y-derivative coefficients, each row is a basis with the monomial coefficients given in the columns
    Eigen::Matrix2Xd _nodes; // Lagrange nodes on a unit right triangle

    Eigen::VectorXd evalPhi (double x, double y) const noexcept; // _phi evaluated at x and y
    Eigen::VectorXd evalPhiX(double x, double y) const noexcept; // _phix evaluated at x and y
    Eigen::VectorXd evalPhiY(double x, double y) const noexcept; // _phiy evaluated at x and y
    Eigen::Matrix2Xd getLagrangeNodes() const noexcept; // compute the Lagrange nodes
};

#endif