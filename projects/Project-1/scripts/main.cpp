// src/main.cpp
#include <iostream>
#include <string>
#include <vector>
#include <filesystem>

#include "../include/ReadGRI.h"
#include "../include/ReadConnData.h"
#include "../include/MeshVerification.h"
#include "TriangularMesh.h"

// declare demo entrypoint (defined elsewhere)
extern int run_refine_demo(int argc, char** argv);

static bool verificationSuite();

int main(int argc, char** argv) {

    // ---------------- CLI dispatcher ----------------
    // Example:
    //   ./aero623 refine <mesh.gri> <bladeupper.txt> <bladelower.txt> <out.gri> [alpha] [maxIters]
    if (argc >= 2) {
        std::string mode = argv[1];

        if (mode == "refine") {
            // Forward all args to refine demo
            return run_refine_demo(argc, argv);
        }

        if (mode == "verify") {
            // fall through to verification below
        } else if (mode != "verify") {
            std::cerr
                << "Unknown mode: " << mode << "\n\n"
                << "Usage:\n"
                << "  " << argv[0] << " verify\n"
                << "  " << argv[0] << " refine <mesh.gri> <bladeupper.txt> <bladelower.txt> <out.gri> [alpha] [maxIters]\n";
            return 2;
        }
    }

    // ---------------- default behavior: verification ----------------
    // If you need this, keep it. Otherwise remove it.
    TriangularMesh mesh("projects/Project-1/mesh_coarse.gri");
    mesh.writeGri("projects/Project-1/mesh_coarse.gri");

    if (!verificationSuite()) {
        throw std::runtime_error("Verification test(s) failed");
    }

    return 0;
}

static bool verificationSuite() {
    std::cout << "Running mesh verification suite..........................\n";

    std::filesystem::path currentDir = std::filesystem::current_path();
    std::cout << "Current working directory: " << currentDir.string() << "\n";

    // ---- Test mesh ----
    {
        std::string gri = "projects/Project-1/test.gri";
        std::vector<std::string> txt = {
            "projects/Project-1/testperiodicEdges.txt",
            "projects/Project-1/testI2E.txt",
            "projects/Project-1/testIn.txt",
            "projects/Project-1/testB2E.txt",
            "projects/Project-1/testBn.txt"
        };
        std::cout << "------------Mesh Verification for Test Grid------------\n";
        meshVerification(gri, txt);
        std::cout << "\n";
    }

    // ---- Coarse mesh ----
    {
        std::string gri = "projects/Project-1/mesh_coarse.gri";
        std::vector<std::string> txt = {
            "projects/Project-1/mesh_coarseperiodicEdges.txt",
            "projects/Project-1/mesh_coarseI2E.txt",
            "projects/Project-1/mesh_coarseIn.txt",
            "projects/Project-1/mesh_coarseB2E.txt",
            "projects/Project-1/mesh_coarseBn.txt"
        };
        std::cout << "------------Mesh Verification for Coarse Grid------------\n";
        meshVerification(gri, txt);
        std::cout << "\n";
    }

    // ---- Local refined mesh ----
    {
        std::string gri = "projects/Project-1/mesh_refined_2394.gri";
        std::vector<std::string> txt = {
            "projects/Project-1/mesh_refined_2394periodicEdges.txt",
            "projects/Project-1/mesh_refined_2394I2E.txt",
            "projects/Project-1/mesh_refined_2394In.txt",
            "projects/Project-1/mesh_refined_2394B2E.txt",
            "projects/Project-1/mesh_refined_2394Bn.txt"
        };
        std::cout << "------------Mesh Verification for Locally Refined Grid------------\n";
        meshVerification(gri, txt);
        std::cout << "\n";
    }

    // ---- Global refined meshes ----
    auto verify_global = [](int k) {
        std::string base = "projects/Project-1/meshGlobalRefined" + std::to_string(k);
        std::string gri  = base + ".gri";
        std::vector<std::string> txt = {
            base + "periodicEdges.txt",
            base + "I2E.txt",
            base + "In.txt",
            base + "B2E.txt",
            base + "Bn.txt"
        };
        std::cout << "------------Mesh Verification for Global Refined Grid " << k << " ------------\n";
        meshVerification(gri, txt);
        std::cout << "\n";
    };

    verify_global(1);
    verify_global(2);
    verify_global(3);

    return true;
}
