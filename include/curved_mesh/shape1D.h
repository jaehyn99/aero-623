#ifndef SHAPEL1D_H
#define SHAPEL1D_H

int shapeL1D(double sig, int q, double **pphi);
int gradientL1D(double sig, int q, double **qgphi);

// Eigen wrappers
Eigen::MatrixXd shapeL1D_quad(const Eigen::VectorXd& xiq, int p);
Eigen::MatrixXd gradL1D_quad(const Eigen::VectorXd& xiq, int p);

#endif