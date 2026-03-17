// // Build the 1D Lagrange points
// // Stand alone executable to build the 1D Lagrange points for curved edge nodes. 
// // This is not meant to be a general-purpose tool,
// // A quick utility to generate the points for our specific blade geometry and mesh.
// // Compile: g++ -std=c++17 src/curved_mesh/geom_snap.cpp src/geom_snap_main.cpp -Iinclude -o geom_snap
// // Run: ./geom_snap

// #include "curved_mesh/geom_snap.h"
// #include <iostream>

// int main() {
//     auto upper = readBladePoints("/home/jaehyn/CFD2/aero-623/projects/Project-3/bladeupper.txt");
//     auto lower = readBladePoints("/home/jaehyn/CFD2/aero-623/projects/Project-3/bladelower.txt");
//     auto lowerShifted = offset(lower, 0.0, 18.0);

//     auto sUpper = arcLength(upper);
//     auto sLower = arcLength(lowerShifted);

//     ParametricSpline2D upperSpline, lowerSpline;
//     upperSpline.sx.build(sUpper, getX(upper));
//     upperSpline.sy.build(sUpper, getY(upper));
//     lowerSpline.sx.build(sLower, getX(lowerShifted));
//     lowerSpline.sy.build(sLower, getY(lowerShifted));

//     GriMesh mesh = readGri("/home/jaehyn/CFD2/aero-623/projects/Project-1/mesh_coarse.gri");

//     int q = 1;  // choose edge order here among 1,2,3

//     auto upperCurved = buildCurvedEdges(
//         "Curve1",
//         mesh.curves["Curve1"],
//         mesh.nodes,
//         upperSpline,
//         sUpper.front(),
//         sUpper.back(),
//         q
//     );

//     auto lowerCurved = buildCurvedEdges(
//         "Curve5",
//         mesh.curves["Curve5"],
//         mesh.nodes,
//         lowerSpline,
//         sLower.front(),
//         sLower.back(),
//         q
//     );

//     writeCurvedEdges("upper_curved_edges.txt", upperCurved);
//     writeCurvedEdges("lower_curved_edges.txt", lowerCurved);

//     return 0;
// }