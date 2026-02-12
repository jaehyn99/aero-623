// #include "SSP_RK2.h"

// StateMesh SSP_RK2::integrate(const Function& f, const StateMesh& u0, double t, double dt) const{
//     StateMesh u1(u0.mesh()), u2(u0.mesh());
//     u1 = u0.array() + dt*f(t, u0).array();
//     u2 = 0.5*u0.array() + 0.5*u1.array() + 0.5*dt*f(t+dt, u1).array();
//     return u2;
// }

// StateMesh SSP_RK2::integrate(const Function& f, const StateMesh& u0, double t, const Eigen::ArrayXd& dt) const{
//     StateMesh u1(u0.mesh()), u2(u0.mesh());
//     u1 = u0.array() + dt*f(t,u0).array();
//     u2 = 0.5*u0.array() + 0.5*u1.array() + 0.5*dt*f(t+dt.minCoeff(), u1).array();
//     return u2;
// }