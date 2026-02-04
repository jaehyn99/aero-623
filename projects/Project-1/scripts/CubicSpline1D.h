#ifndef CUBIC_SPLINE_1D_H
#define CUBIC_SPLINE_1D_H

#include <Eigen/Dense>

namespace mesh {

// Stores knots Si, values Xi, and derivatives dXi for Hermite-like cubic spline.
class CubicSpline1D {
public:
    CubicSpline1D() = default;

    // Set spline from knots/values/derivatives
    void set(const Eigen::VectorXd& S,
             const Eigen::VectorXd& X,
             const Eigen::VectorXd& dX);

    // Evaluate x(s) (Python: splineFun)
    double eval(double s) const;

    // Evaluate dx/ds (Python: diffSplineFun)  [Optional but implemented]
    double deriv(double s) const;

    const Eigen::VectorXd& knots() const { return S_; }
    const Eigen::VectorXd& values() const { return X_; }
    const Eigen::VectorXd& derivs() const { return dX_; }

private:
    Eigen::VectorXd S_, X_, dX_;

    // Find interval index i such that S[i] <= s <= S[i+1]
    Eigen::Index findInterval(double s) const;
};

} // namespace mesh

#endif
