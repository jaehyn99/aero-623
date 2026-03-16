#pragma once
#include <Eigen/Dense>
#include "shapeFunctions/shapeFunctions.hpp"
#include "quadrature/femQuadrature.hpp"
#include <iostream>

class referenceElement {
	public:
		Eigen::MatrixXd phi;
		Eigen::MatrixXd phiXi;
		Eigen::MatrixXd phiEta;
		Eigen::MatrixXd xiQ;
		Eigen::MatrixXd xiL;
		Eigen::MatrixXd MRef;
		Eigen::VectorXd wQ;
		int p, q, nQ, nL;
		referenceElement(int p_, int q_) : p(p_), q(q_) {
			p = p_; q = q_;
			shapeFunctions shape(p_);
			femQuadrature quad(q_);
			xiL = shape.getRefLagrangePoints(p_);
			Eigen::MatrixXd phiCoeffs = shape.getShapeFuncCoeffs(p_); // Gets shape function and derivative coefficients
			Eigen::MatrixXd phiXiCoeffs = shape.getShapeFuncXiCoeffs(p_);
			Eigen::MatrixXd phiEtaCoeffs = shape.getShapeFuncEtaCoeffs(p_);

			xiQ = quad.getQuadXi(q_); // Gets quadrature points in reference space and weights
			wQ = quad.getQuadW(q_);
			nQ = wQ.size();
			nL = (p_ + 1)*(p_ + 2)/2; // Computes number of Lagrange points per triangle
			Eigen::VectorXd mon(nL); 
			phi.resize(nQ, nL);
			phiXi.resize(nQ, nL);
			phiEta.resize(nQ, nL);
			MRef.resize(nL, nL);
			
			phi.setZero(); phiXi.setZero(); phiEta.setZero(); MRef.setZero();

			// Computes monomial coefficients
			for (int ii = 0; ii < nQ; ii++) {
				for (int jj = 0; jj < p_ + 1; jj++) {
					for (int kk = 0; kk < jj + 1; kk++) {
						mon(jj*(jj + 1)/2 + kk) = std::pow(xiQ(ii, 0), jj - kk)*std::pow(xiQ(ii, 1), kk);
					}
				}

				// computes shape functions and derivatives at quadrature points
				for (int jj = 0; jj < nL; jj++) {
					for (int kk = 0; kk < nL; kk++) {
						phi(ii, jj) = phi(ii, jj) + phiCoeffs(jj, kk)*mon(kk);
					}
					for (int kk = 0; kk < p_*(p_ + 1)/2; kk++) {
						phiXi(ii, jj) = phiXi(ii, jj) + phiXiCoeffs(kk, jj)*mon(kk);
						phiEta(ii, jj) = phiEta(ii, jj) + phiEtaCoeffs(kk, jj)*mon(kk);
					}
				}
			}

			// Computes reference mass matrix using quadrature
			for (int ii = 0; ii < nL; ii++) {
				for (int jj = 0; jj < nL; jj++) {
					for (int kk = 0; kk < nQ; kk++) {
						MRef(ii, jj) = MRef(ii, jj) + wQ(kk)*phi(kk, ii)*phi(kk, jj);
					}
				}
			}
		}
	protected:
		int p_;
		int q_;
};
