#pragma once

#include <array>
#include <cstddef>
#include <memory>
#include <string>
#include <vector>


#include "mesh/TriangularMesh.h"

class FirstorderEuler {
public:
    struct SolverConfig {
        double gamma = 1.4;
        double gasConstant = 1.0;

        double rho0 = 1.0;
        double a0 = 1.0;
        // Stagnation pressure is computed as p0 = rho0 * a0^2 / gamma.
        double p0 = 1.0;
        double alpha = 0.0;
        double pout = 1.0;

        double cfl = 0.1;
        double finalTime = 2.0;
        std::size_t maxIterations = 100000;
        bool localTimeStepping = true;

        double residualTolerance = 1e-5;
        std::size_t saveEvery = 1000;
        std::string outputPrefix = "sol";
        std::string fluxScheme = "roe"; // "roe" or "hlle"

        double initialMach = 0.0;

        // Boundary-curve role mapping (periodic remains 2<->4 and 6<->8).
        int inflowCurve = 3;
        int outflowCurve = 7;

        // Optional debug guard: validate mesh/connectivity arrays after load.
        bool validateMeshOnLoad = true;

        // Lightweight runtime diagnostics (printed every debugEvery iterations).
        bool enableDebugPrints = true;
        std::size_t debugEvery = 10;

        // One-line progress table cadence ("it t dt ||R||2").
        std::size_t statusEvery = 1000;
    };

    struct MeshInputs {
        std::string meshFile;
    };

    using Conserved = std::array<double, 4>; // [rho, rhou, rhov, rhoE]
    using Vec2 = std::array<double, 2>;
    Vec2 cellCentroid(std::size_t ei) const;


    struct InteriorFace {
        std::size_t elemL = 0;
        std::size_t faceL = 0;
        std::size_t elemR = 0;
        std::size_t faceR = 0;
        Vec2 normal{0.0, 0.0};
        Vec2 center{0.0, 0.0};
        double length = 0.0;
    };

    struct BoundaryFace {
        std::size_t elem = 0;
        std::size_t localFace = 0;
        std::size_t boundaryGroup = 0;
        std::string boundaryTitle;
        Vec2 normal{0.0, 0.0};
        Vec2 center{0.0, 0.0};
        double length = 0.0;
    };

    struct EdgeFluxContribution {
        std::size_t ownerElem = 0;
        std::size_t neighborElem = 0;
        Vec2 normal{0.0, 0.0};
        double edgeLength = 0.0;
        Conserved flux{0.0, 0.0, 0.0, 0.0};
        double spectralRadius = 0.0;
    };

    enum class BoundaryKind {
        InflowSteady,
        InflowUnsteady,
        OutflowSubsonic,
        WallSlip,
        Periodic
    };

    FirstorderEuler(MeshInputs inputs, SolverConfig config);

    // Part 1 API: mesh processing + conservative initialization.
    void loadInputs();
    // Build one uniform initial conservative state from inlet stagnation data
    // and the initial Mach guess, then copy it to every cell.
    void initUniformState();

    
    // Later parts (kept but currently skeletons).
    void advanceToConvergedOrFinalTime();

    // Run modes for first-order solver collaboration:
    // - steady with global dt
    // - steady with local dt
    // - unsteady with global dt
    void runSteadyGlobal();
    void runSteadyLocal();
    void runUnsteadyGlobal();

    // Export solver state as legacy VTK (cell fields include pressure, Mach and Cp).
    void writeSolutionVtk(const std::string& filePath) const;

    const std::vector<Conserved>& solution() const { return U_; }
    const std::vector<Conserved>& residual() const { return residual_; }

private:
    // Part 1: read mesh and build solver-facing arrays.
    void readMeshAndConnectivity();
    void buildFacesFromTriangularMesh();
    void computePerimeterFromMesh();

    void validateLoadedArrays() const; // debug-only, gated by config_.validateMeshOnLoad

    // Part 2+ placeholders.
    std::vector<EdgeFluxContribution> computeEdgeFluxesAndWaveSpeeds() const;
    Conserved computeBoundaryFluxFromModules(const BoundaryFace& f,
                                              const Conserved& UL,
                                              BoundaryKind kind) const;
    void assembleResidualFromEdgeFluxes(const std::vector<EdgeFluxContribution>& edges);
    double computeGlobalDtFromWaveSpeeds(const std::vector<EdgeFluxContribution>& edges) const;
    std::vector<double> computeLocalDtFromWaveSpeeds(const std::vector<EdgeFluxContribution>& edges) const;
    void updateStateGlobalDt(double dt);
    void updateStateLocalDt(const std::vector<double>& dtLocal);
    void resetMarchState();
    void advance(bool stopByTime);
    static double l2Norm(const std::vector<Conserved>& values);
    double cellPressure(const Conserved& U) const;
    void printBoundaryConditionSummary() const;
    void printIterationDiagnostics(const std::vector<EdgeFluxContribution>& edges,
                                   const std::vector<double>* dtLocal,
                                   double dtUsed,
                                   double normR) const;
    double maxWallBoundaryMassFlux() const;

    MeshInputs inputs_;
    SolverConfig config_;
    std::unique_ptr<TriangularMesh> triMesh_;

    // Mesh data
    std::vector<Vec2> nodes_;
    std::vector<std::array<std::size_t, 3>> elements_;
    std::vector<InteriorFace> interiorFaces_;
    std::vector<BoundaryFace> boundaryFaces_;
    std::vector<std::pair<std::size_t, std::size_t>> periodicEdges_;

    // Geometry scalars per element.
    std::vector<double> area_;
    std::vector<double> perimeter_;

    // Conservative state per element.
    std::vector<Conserved> U_;
    std::vector<Conserved> residual_;

    double time_ = 0.0;
    std::size_t iteration_ = 0;
};
