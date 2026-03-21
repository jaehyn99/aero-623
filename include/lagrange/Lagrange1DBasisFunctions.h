// The set of 2D Lagrange basis functions defined on a unit right triangle

#ifndef LAGRANGE_1D_BASIS_FUNCTIONS_H
#define LAGRANGE_1D_BASIS_FUNCTIONS_H

#include "Eigen/Dense"
class Lagrange1DBasisFunctions{
    public:
    Lagrange1DBasisFunctions(int p);

    int p() const noexcept{ return _p; }
    int Np() const noexcept{ return _p+1; } 
    Eigen::MatrixXd phi() const noexcept{ return _phi; }
    Eigen::MatrixXd phix() const noexcept{ return _phix; }
    Eigen::VectorXd nodes() const noexcept{ return Eigen::VectorXd::LinSpaced(_p+1, 0.0, 1.0); }

    // Evaluation of the basis functions and their derivatives
    Eigen::VectorXd evalPhi (double x) const noexcept; // _phi evaluated at x and y
    Eigen::VectorXd evalPhiX(double x) const noexcept; // _phix evaluated at x and y

    // Evaluation of a function (given by the associated weights) and its derivatives
    double funcEval (double x, const Eigen::VectorXd& coeff);
    double funcXEval(double x, const Eigen::VectorXd& coeff);

    // Multi-D overloads that returns an Ns-by-1 vector
    // coeff is of dimension Ns*Np, where Ns = number of states, Np = number of basis functions;
    Eigen::VectorXd funcEval (double x, const Eigen::MatrixXd& coeff);
    Eigen::VectorXd funcXEval(double x, const Eigen::MatrixXd& coeff);

    protected:
    int _p; // polynomial order
    Eigen::MatrixXd _phi; // matrix of basis function coefficients, each row is a basis with the monomial coefficients given in the columns
    Eigen::MatrixXd _phix; // matrix of derivative coefficients, each row is a basis with the monomial coefficients given in the columns
};

#endif