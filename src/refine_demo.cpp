// src/refine_demo.cpp
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
#include <Eigen/Dense>

#include "mesh/BladeGeometry.h"
#include "mesh/Projection2D.h"
#include "mesh/SizingFunction.h"
#include "mesh/LocalRefinement.h"
#include "mesh/GriIO.h"

static std::string replaceExt(const std::string& path, const std::string& newExt)
{
    const auto pos = path.find_last_of('.');
    if (pos == std::string::npos) return path + newExt;
    return path.substr(0, pos) + newExt;
}

int run_refine_demo(int argc, char** argv)
{
    // Usage:
    //   main.exe refine <mesh.gri> <bladeupper.txt> <bladelower.txt> <out.gri> [alpha] [maxIters]
    if (argc < 6) {
        std::cerr << "Usage: " << argv[0]
                  << " refine <mesh.gri> <bladeupper.txt> <bladelower.txt> <out.gri> [alpha] [maxIters]\n";
        return 2;
    }

    const std::string meshFile  = argv[2];
    const std::string upperFile = argv[3];
    const std::string lowerFile = argv[4];
    const std::string outFile   = argv[5];

    const double alpha   = (argc >= 7) ? std::stod(argv[6]) : 1.2;
    const int    maxIters = (argc >= 8) ? std::stoi(argv[7]) : 5;

    mesh::GriMesh2D meshIn = mesh::readGri2D(meshFile);

    mesh::BladeGeometry blade;
    auto mesh = mesh::readGri2D(meshFile);
    double Ly = mesh.V.col(1).maxCoeff() - mesh.V.col(1).minCoeff();

    // if the lower blade should align to the “upper wall” copy:
    blade.loadAndFitTwo(upperFile, lowerFile, 1e-12, 18);

    mesh::Projection2D projector;

    // Sizing parameters
    const double hMin = 0.02;
    const double hMax = 1.20; // at least 0.9
    const double d0   = 1.5;

    // Compute LE/TE robustly from blade sampling
    auto bladeXMinMax = [&](int Ns=400) {
        double xmin =  1e300, xmax = -1e300;
        for (int i=0; i<Ns; ++i) {
            const double tu = double(i) / double(Ns-1);
            const double tl = tu;
            const double su = blade.sUmin() + (blade.sUmax()-blade.sUmin()) * tu;
            const double sl = blade.sLmin() + (blade.sLmax()-blade.sLmin()) * tl;
            const auto pu = blade.evalUpper(su);
            const auto pl = blade.evalLower(sl);
            xmin = std::min(xmin, std::min(pu.x(), pl.x()));
            xmax = std::max(xmax, std::max(pu.x(), pl.x()));
        }
        return std::pair<double,double>(xmin, xmax);
    };

    const auto [xLE, xTE] = bladeXMinMax();
    const double xSig = 0.1;

    mesh::SizingFunction sizeFun(hMin, hMax, d0, xLE, xTE, xSig, 0.6, 2.0);

    // --- refine ---
    mesh::GriMesh2D meshOut = mesh::LocalRefinement::refineMesh(
        meshIn, blade, projector, sizeFun, alpha, maxIters
    );

    auto testFar = [&](double x, double y){
    Eigen::Vector2d p(x,y);
    auto r = mesh::LocalRefinement::projectDistance(p, meshIn.V, blade, projector);
    std::cout << "p=("<<x<<","<<y<<") dist="<<r.dist<<" curveId="<<r.curveId
              <<" xProj=("<<r.xProj.x()<<","<<r.xProj.y()<<")\n";
    };

    double xmin = meshIn.V.col(0).minCoeff(), xmax = meshIn.V.col(0).maxCoeff();
    double ymin = meshIn.V.col(1).minCoeff(), ymax = meshIn.V.col(1).maxCoeff();

    testFar(xmin, ymin);
    testFar(xmin, ymax);
    testFar(xmax, ymin);
    testFar(xmax, ymax);
    // Debug

    mesh::writeGri2D(outFile, meshOut);

    // Compute fields on refined mesh
    Eigen::VectorXd distNode, hNode;
    mesh::LocalRefinement::computeNodeDistanceAndSize(meshOut.V, blade, projector, sizeFun, distNode, hNode);

    const std::string distPath = replaceExt(outFile, ".walldist.txt");
    const std::string hPath    = replaceExt(outFile, ".hnode.txt");

    {
        std::ofstream dfile(distPath);
        if (!dfile) throw std::runtime_error("cannot write: " + distPath);
        dfile.setf(std::ios::scientific);
        dfile.precision(17);
        for (int i = 0; i < distNode.size(); ++i) dfile << distNode(i) << "\n";
    }
    {
        std::ofstream hfile(hPath);
        if (!hfile) throw std::runtime_error("cannot write: " + hPath);
        hfile.setf(std::ios::scientific);
        hfile.precision(17);
        for (int i = 0; i < hNode.size(); ++i) hfile << hNode(i) << "\n";
    }

    std::cout << "Refinement done.\n"
              << "Input : nodes=" << meshIn.V.rows()  << " elems=" << meshIn.E.rows()  << "\n"
              << "Output: nodes=" << meshOut.V.rows() << " elems=" << meshOut.E.rows() << "\n"
              << "Wrote: " << outFile << "\n"
              << "Wrote: " << distPath << "\n"
              << "Wrote: " << hPath << "\n";

    return 0;
}
