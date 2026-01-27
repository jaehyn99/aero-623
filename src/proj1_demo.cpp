// src/proj1_demo.cpp
#include <iostream>
#include <iomanip>
#include <string>
#include <Eigen/Dense>

#include "mesh/BladeGeometry.h"
#include "mesh/Projection2D.h"

int run_proj1_demo(int argc, char** argv) {
    // Usage:
    //   main.exe proj1 <bladeupper.txt> <bladelower.txt> [px py]
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0]
                  << " proj1 <bladeupper.txt> <bladelower.txt> [px py]\n";
        return 2;
    }

    const std::string upper = argv[2];
    const std::string lower = argv[3];

    double px = 15.0, py = -5.0;
    if (argc >= 6) {
        px = std::stod(argv[4]);
        py = std::stod(argv[5]);
    }

    mesh::BladeGeometry blade;
    blade.loadAndFitTwo(upper, lower);

    mesh::Projection2D projector;

    Eigen::Vector2d p(px, py);
    auto res = projector.projectToBlade(blade, p);

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "p = [" << p.x() << ", " << p.y() << "]\n";
    std::cout << "s* = " << res.s << "\n";
    std::cout << "proj = [" << res.xProj.x() << ", " << res.xProj.y() << "]\n";
    std::cout << "dist = " << res.dist << "\n";

    return 0;
}
