#pragma once
#include <Eigen/Dense>
#include "mesh/referenceElement.hpp"
#include "ShapeFunctions/ShapeFunctions.hpp"
#include "quadrature/femQuadrature.hpp"
#include <iostream>

class globalLinearElement {
	public:
		Eigen::MatrixXd xL, M, J, invJ;
		double detJ;
		globalLinearElement(Eigen::MatrixXd& V, referenceElement& refElem) {
			int nL = refElem.nL; int nQ = refElem.nQ;
			Eigen::MatrixXd xiL = refElem.xiL;
			xL.resize(nL, 2);
			J.resize(2, 2);

			// Jacobian for geometrically linear elements
			J << V(1, 0) - V(0, 0), V(2, 0) - V(0, 0),
			     V(1, 1) - V(0, 1), V(2, 1) - V(0, 1);

			// Each global linear element stores a unique Jacobian based on vertex locations and 
			// mass matrix (reference element's mass matrix scaled by Jacobian determinant)
			detJ = J.determinant();
			invJ = J.inverse();
			M = refElem.MRef*detJ;

			// Also Lagrange nodes in global space (not correctly computed yet)
/*
			for (int ii = 0; ii < nL; ii++) {
				for (int jj = 0; jj < nL; jj++) {
			        	xL(ii, 0) = V(0, jj) + refElem.xiL(ii,0)*(V(1, jj) - V(0, jj)) + refElem.xiL(ii,1)*(V(2, jj) - V(0, jj));
			        	xL(ii, 1) = V(0, jj) + refElem.xiL(ii,0)*(V(1, jj) - V(0, jj)) + refElem.xiL(ii,1)*(V(2, jj) - V(0, jj));
				}
			}
*/		}
};

