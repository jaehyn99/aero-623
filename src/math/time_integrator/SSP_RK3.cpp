// #include "SSP_RK3.h"

// StateMesh SSP_RK3::integrate(const Function& f, const StateMesh& u0, double t, double dt) const{
//     StateMesh u1(u0.mesh()), u2(u0.mesh()), u3(u0.mesh());
//     u1 = u0.array() + dt*f(t, u0).array();
//     u2 = 0.75*u0.array() + 0.25*u1.array() + 0.25*dt*f(t+dt, u1).array();
//     u3 = u0.array()/3 + u2.array()*(2.0/3) + 2.0/3*dt*f(t+0.5*dt, u2).array();
//     return u3;
// }

// StateMesh SSP_RK3::integrate(const Function& f, const StateMesh& u0, double t, const Eigen::ArrayXd& dt) const{
//     StateMesh u1(u0.mesh()), u2(u0.mesh()), u3(u0.mesh());
//     u1 = u0.array() + dt*f(t, u0).array();
//     u2 = 0.75*u0.array() + 0.25*u1.array() + 0.25*dt*f(t+dt.minCoeff(), u1).array();
//     u3 = u0.array()/3 + u2.array()*(2.0/3) + 2.0/3*dt*f(t+0.5*dt.minCoeff(), u2).array();
//     return u3;
// }