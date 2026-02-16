#include "solver/BoundaryConditions.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace {

int parseCurveId(const std::string& title) {
    std::string digits;
    for (char c : title) {
        if (std::isdigit(static_cast<unsigned char>(c))) {
            digits.push_back(c);
        }
    }
    if (digits.empty()) {
        return -1;
    }
    return std::stoi(digits);
}

} // namespace

EulerBoundaryConditions::EulerBoundaryConditions(Config config)
    : config_(std::move(config)) {}

EulerBoundaryConditions::Type EulerBoundaryConditions::typeFromCurveTitle(const std::string& title) const {
    const int curveID = parseCurveId(title);

    // Project convention from mesh titles.
    if (curveID == 2 || curveID == 4 || curveID == 6 || curveID == 8) {
        return Type::Periodic;
    }
    if (curveID == 1 || curveID == 5) {
        return Type::WallSlip;
    }
    if (curveID == config_.outflowCurve) {
        return Type::OutflowSubsonic;
    }
    if (curveID == config_.inflowCurve) {
        return Type::InflowSteady;
    }

    // Fallback for historical meshes using Curve3/Curve7 pair.
    if (curveID == 3 || curveID == 7) {
        return (curveID == 3) ? Type::InflowSteady : Type::OutflowSubsonic;
    }
    return Type::InflowSteady;
}

EulerBoundaryConditions::Conserved EulerBoundaryConditions::boundaryState(
    Type type,
    const Conserved& UI,
    const Vec2& n,
    const Context& ctx,
    const Conserved* periodicPartner) const {
    switch (type) {
    case Type::InflowSteady:
        return inflowSteadyState(UI, n, config_.pt);
    case Type::InflowUnsteady:
        return inflowUnsteadyState(UI, n, ctx);
    case Type::OutflowSubsonic:
        return outflowSubsonicState(UI, n);
    case Type::WallSlip:
        return wallSlipState(UI, n);
    case Type::Periodic:
        if (periodicPartner == nullptr) {
            throw std::runtime_error("Periodic boundary requires partner state.");
        }
        return *periodicPartner;
    }

    return UI;
}

EulerBoundaryConditions::Vec2 EulerBoundaryConditions::normalize(Vec2 n) {
    const double norm = std::sqrt(n[0] * n[0] + n[1] * n[1]);
    if (norm <= 1e-14) {
        throw std::runtime_error("Boundary normal has near-zero magnitude.");
    }
    n[0] /= norm;
    n[1] /= norm;
    return n;
}

double EulerBoundaryConditions::dot(const Vec2& a, const Vec2& b) {
    return a[0] * b[0] + a[1] * b[1];
}

double EulerBoundaryConditions::pressure(const Conserved& U) const {
    const double rho = std::max(1e-14, U[0]);
    const double u = U[1] / rho;
    const double v = U[2] / rho;
    return (config_.gamma - 1.0) * (U[3] - 0.5 * rho * (u * u + v * v));
}

EulerBoundaryConditions::Conserved EulerBoundaryConditions::wallSlipState(const Conserved& UI, const Vec2& nRaw) const {
    const Vec2 n = normalize(nRaw);

    Conserved Ub = UI;
    const double rho = std::max(1e-14, UI[0]);
    const double u = UI[1] / rho;
    const double v = UI[2] / rho;
    const double p = std::max(1e-14, pressure(UI));

    const double un = u * n[0] + v * n[1];
    const double uw = u - 2.0 * un * n[0];
    const double vw = v - 2.0 * un * n[1];

    Ub[1] = rho * uw;
    Ub[2] = rho * vw;
    Ub[3] = p / (config_.gamma - 1.0) + 0.5 * rho * (uw * uw + vw * vw);
    return Ub;
}

EulerBoundaryConditions::Conserved EulerBoundaryConditions::inflowSteadyState(
    const Conserved& UI,
    const Vec2& nRaw,
    double ptOverride) const {
    const Vec2 n = normalize(nRaw);

    const double gam = config_.gamma;
    const double R = config_.gasConstant;
    const Vec2 nin{std::cos(config_.alpha), std::sin(config_.alpha)};
    const double dn = dot(nin, n);

    const double rI = std::max(1e-14, UI[0]);
    const Vec2 VI{UI[1] / rI, UI[2] / rI};
    const double unI = dot(VI, n);
    const double pI = std::max(1e-14, (gam - 1.0) * (UI[3] - 0.5 * rI * dot(VI, VI)));
    const double cI = std::sqrt(gam * pI / rI);

    const double Jp = unI + 2.0 * cI / (gam - 1.0);

    const double A = gam * R * config_.Tt * dn * dn - 0.5 * (gam - 1.0) * Jp * Jp;
    const double B = 4.0 * gam * R * config_.Tt * dn / (gam - 1.0);
    const double C = 4.0 * gam * R * config_.Tt / ((gam - 1.0) * (gam - 1.0)) - Jp * Jp;

    double Mb = 0.1;
    if (std::abs(A) > 1e-14) {
        const double disc = B * B - 4.0 * A * C;
        if (disc >= 0.0) {
            const double sqrtDisc = std::sqrt(disc);
            const double m1 = (-B - sqrtDisc) / (2.0 * A);
            const double m2 = (-B + sqrtDisc) / (2.0 * A);
            const bool m1Pos = m1 > 0.0;
            const bool m2Pos = m2 > 0.0;
            if (m1Pos && m2Pos) {
                Mb = std::min(m1, m2);
            } else if (m1Pos) {
                Mb = m1;
            } else if (m2Pos) {
                Mb = m2;
            }
        }
    }

    const double Tb = config_.Tt / (1.0 + 0.5 * (gam - 1.0) * Mb * Mb);
    const double pb = ptOverride * std::pow(Tb / config_.Tt, gam / (gam - 1.0));
    const double rb = std::max(1e-14, pb / (R * Tb));
    const double cb = std::sqrt(gam * pb / rb);
    const Vec2 Vb{Mb * cb * nin[0], Mb * cb * nin[1]};
    const double rEb = pb / (gam - 1.0) + 0.5 * rb * dot(Vb, Vb);

    Conserved Ub{};
    Ub[0] = rb;
    Ub[1] = rb * Vb[0];
    Ub[2] = rb * Vb[1];
    Ub[3] = rEb;
    return Ub;
}

EulerBoundaryConditions::Conserved EulerBoundaryConditions::inflowUnsteadyState(
    const Conserved& UI,
    const Vec2& n,
    const Context& ctx) const {
    const double yRot = ctx.faceCenter[1] + config_.yRot;
    const double yStator = yRot + config_.Vrot * ctx.time;
    const double z = yStator / std::max(1e-14, config_.deltaY);
    const double eta = z - std::floor(z) - 0.5;

    const double shape = std::exp(-(eta * eta) / (2.0 * config_.delta * config_.delta));
    const double rho0Local = config_.rho0 * (1.0 - config_.fwake * shape);

    // Keep Tt fixed and adjust pt with stagnation density relation p_t = rho_t * R * T_t.
    const double ptLocal = rho0Local * config_.gasConstant * config_.Tt;
    return inflowSteadyState(UI, n, ptLocal);
}

EulerBoundaryConditions::Conserved EulerBoundaryConditions::outflowSubsonicState(
    const Conserved& UI,
    const Vec2& nRaw) const {

    const Vec2 n = normalize(nRaw);
    const double gam = config_.gamma;

    const double rhoI = std::max(1e-14, UI[0]);
    const Vec2 VI{UI[1] / rhoI, UI[2] / rhoI};
    const double unI = dot(VI, n);
    const double pI  = std::max(1e-14, pressure(UI));
    const double cI  = std::sqrt(gam * pI / rhoI);

    // ---- NEW: handle regime robustness ----
    // 1) Supersonic outflow: DO NOT impose pressure; extrapolate
    if (unI >= cI) {
        return UI;
    }

    // 2) Backflow at outflow boundary: must prescribe inflow-like data
    // A simple robust choice: use farfield inflow state (steady) or clamp normal velocity to 0.
    if (unI < 0.0) {
        // Option A (recommended): treat as inflow with total conditions
        return inflowSteadyState(UI, n, config_.pt);

        // Option B (less physical but sometimes stabilizes): prevent inflow
        // Conserved Ub = UI;
        // const Vec2 Vt{VI[0] - unI*n[0], VI[1] - unI*n[1]};
        // const double unb = 0.0;
        // const Vec2 Vb{Vt[0] + unb*n[0], Vt[1] + unb*n[1]};
        // Ub[1] = rhoI*Vb[0];
        // Ub[2] = rhoI*Vb[1];
        // Ub[3] = pI/(gam-1.0) + 0.5*rhoI*dot(Vb,Vb);
        // return Ub;
    }

    // ---- Existing subsonic outflow logic (pressure specified) ----
    const double Jp = unI + 2.0 * cI / (gam - 1.0);
    const double Splus = pI / std::pow(rhoI, gam);

    const double pb = std::max(1e-14, config_.pout);
    const double rhob = std::max(1e-14, std::pow(pb / Splus, 1.0 / gam));
    const double cb = std::sqrt(gam * pb / rhob);
    const double unb = Jp - 2.0 * cb / (gam - 1.0);

    const Vec2 VbTangent{VI[0] - unI * n[0], VI[1] - unI * n[1]};
    const Vec2 Vb{VbTangent[0] + unb * n[0], VbTangent[1] + unb * n[1]};

    Conserved Ub{};
    Ub[0] = rhob;
    Ub[1] = rhob * Vb[0];
    Ub[2] = rhob * Vb[1];
    Ub[3] = pb / (gam - 1.0) + 0.5 * rhob * dot(Vb, Vb);
    return Ub;
}
