#pragma once

#include <array>
#include <cstddef>
#include <string>

class EulerBoundaryConditions {
public:
    using Conserved = std::array<double, 4>; // [rho, rhou, rhov, rhoE]
    using Vec2 = std::array<double, 2>;

    enum class Type {
        InflowSteady,
        InflowUnsteady,
        OutflowSubsonic,
        WallSlip,
        Periodic
    };

    struct Config {
        double gamma = 1.4;
        double gasConstant = 1.0;

        // Inflow total conditions.
        double Tt = 1.0;
        double pt = 1.0;
        double alpha = 0.0;

        // Unsteady inflow wake model.
        double rho0 = 1.0;
        double fwake = 0.1;
        double delta = 0.1;
        double deltaY = 18e-3;
        double yRot = 0.0;
        double Vrot = 1.0;

        // Outflow target pressure.
        double pout = 1.0;
    };

    struct Context {
        double time = 0.0;
        Vec2 faceCenter{0.0, 0.0};
    };

    explicit EulerBoundaryConditions(Config config);

    Conserved boundaryState(Type type,
                            const Conserved& UI,
                            const Vec2& n,
                            const Context& ctx,
                            const Conserved* periodicPartner = nullptr) const;

    static Type typeFromCurveTitle(const std::string& title);

private:
    Config config_;

    static Vec2 normalize(Vec2 n);
    static double dot(const Vec2& a, const Vec2& b);

    double pressure(const Conserved& U) const;
    Conserved wallSlipState(const Conserved& UI, const Vec2& n) const;
    Conserved inflowSteadyState(const Conserved& UI, const Vec2& n, double ptOverride) const;
    Conserved inflowUnsteadyState(const Conserved& UI, const Vec2& n, const Context& ctx) const;
    Conserved outflowSubsonicState(const Conserved& UI, const Vec2& n) const;
};
