#include <iostream>
#include <stdexcept>
#include <string>

#include "solver/FirstorderEuler.h"

namespace {

void printUsage(const char* exe) {
    std::cout << "Usage:\n"
              << "  " << exe << " --mode <steady-global|steady-local|unsteady-global> [--mesh <mesh.gri>]\n"
              << "  " << exe << " --mode <steady-global|steady-local|unsteady-global> [--prefix <mesh_prefix>]\n"
              << "\nExamples:\n"
              << "  " << exe << " --mode steady-local --prefix projects/Project-1/mesh_refined_2394\n"
              << "  " << exe << " --mode steady-global --mesh projects/Project-1/mesh_coarse.gri\n";
}

} // namespace

int main(int argc, char** argv) {
    std::string mode = "steady-local";
    std::string meshFile;
    std::string meshPrefix = "projects/Project-1/mesh_refined_2394";

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--mode" && i + 1 < argc) {
            mode = argv[++i];
        } else if (arg == "--mesh" && i + 1 < argc) {
            meshFile = argv[++i];
        } else if (arg == "--prefix" && i + 1 < argc) {
            meshPrefix = argv[++i];
        } else if (arg == "--help" || arg == "-h") {
            printUsage(argv[0]);
            return 0;
        } else {
            std::cerr << "Unknown argument: " << arg << "\n";
            printUsage(argv[0]);
            return 1;
        }
    }

    try {
        if (meshFile.empty()) {
            meshFile = meshPrefix + ".gri";
        }

        FirstorderEuler::MeshInputs inputs;
        inputs.meshFile = meshFile;

        FirstorderEuler::SolverConfig cfg;
        if (mode == "steady-global") {
            cfg.localTimeStepping = false;
            cfg.cfl = 0.1;
            cfg.finalTime = 1e12; // Steady pseudo-time: stop by residual/maxIterations.
        } else if (mode == "steady-local") {
            cfg.localTimeStepping = true;
            cfg.cfl = 0.02;
            cfg.finalTime = 1e12; // Steady pseudo-time: stop by residual/maxIterations.
        } else if (mode == "unsteady-global") {
            cfg.localTimeStepping = false;
            cfg.cfl = 0.1;
            cfg.finalTime = 2.0; // Physical final time for unsteady run.
        } else {
            throw std::runtime_error("Unsupported mode: " + mode);
        }

        // Use legacy HLLE implementation from hlleFlux.hpp for this debug campaign.
        cfg.fluxScheme = "hlle";
        // Write VTK snapshots every 10 iterations to localize discontinuity growth.
        cfg.saveEvery = 10;

        FirstorderEuler solver(inputs, cfg);
        solver.loadInputs();
        solver.initUniformState();

        if (mode == "steady-global") {
            solver.runSteadyGlobal();
        } else if (mode == "steady-local") {
            solver.runSteadyLocal();
        } else {
            solver.runUnsteadyGlobal();
        }
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
