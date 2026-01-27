// src/main.cpp

#include <iostream>
#include <string>

extern int run_proj1_demo(int argc, char** argv);
extern int run_refine_demo(int argc, char** argv);

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr
          << "Usage:\n"
          << "  " << argv[0] << " proj1  <bladeupper.txt> <bladelower.txt> [px py]\n"
          << "  " << argv[0] << " refine <mesh.gri> <bladeupper.txt> <bladelower.txt> <out.gri> [alpha] [maxIters]\n";
        return 2;
    }

    const std::string mode = argv[1];
    if (mode == "proj1")  return run_proj1_demo(argc, argv);
    if (mode == "refine") return run_refine_demo(argc, argv);

    std::cerr << "Unknown mode: " << mode << "\n";
    return 2;
}