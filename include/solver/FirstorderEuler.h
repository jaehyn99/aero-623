#pragma once

#include <array>
#include <cstddef>
#include <memory>
#include <string>
#include <vector>


#include "mesh/TriangularMesh.h"
#include "solver/BoundaryConditions.h"

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

        double cfl = 0.8;
        double finalTime = 2.0;
        std::size_t maxIterations = 20000;
        bool localTimeStepping = true;

        double residualTolerance = 1e-8;
        std::size_t saveEvery = 0;
        std::string outputPrefix = "sol";
        std::string fluxScheme = "roe"; // "roe" or "hlle"

        double initialMach = 0.1;

        // Optional debug guard: validate mesh/connectivity arrays after load.
        bool validateMeshOnLoad = false;
    };

    struct MeshInputs {
        std::string meshFile;
    };

    using Conserved = std::array<double, 4>; // [rho, rhou, rhov, rhoE]
    using Vec2 = std::array<double, 2>;

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
    void assembleResidualFromEdgeFluxes(const std::vector<EdgeFluxContribution>& edges);
    double computeGlobalDtFromWaveSpeeds(const std::vector<EdgeFluxContribution>& edges) const;
    std::vector<double> computeLocalDtFromWaveSpeeds(const std::vector<EdgeFluxContribution>& edges) const;
    void updateStateGlobalDt(double dt);
    void updateStateLocalDt(const std::vector<double>& dtLocal);
    void resetMarchState();
    void advance(bool stopByTime);
    static double l2Norm(const std::vector<Conserved>& values);

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
