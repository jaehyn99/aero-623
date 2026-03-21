#ifndef GEOM_SNAP_H
#define GEOM_SNAP_H

#include <vector>
#include <string>
#include <map>
#include "curved_mesh/geom_snap.h"

struct Point2D {
    double x, y;
};

struct BoundaryPair {
    int n1, n2;   // 1-based node ids
};

struct GriMesh {
    std::vector<Point2D> nodes;   // 1-based in file, stored 0-based here
    std::map<std::string, std::vector<BoundaryPair>> curves;
};

struct CurvedBoundaryEdge {
    std::string curveName;
    int n1, n2;
    std::vector<Point2D> qnodes;
};

// Read blade geometry points
std::vector<Point2D> readBladePoints(const std::string& filename);

// Offset coordinates
std::vector<Point2D> offset(const std::vector<Point2D>& pts, double dx, double dy);

// Arc-length parameterization
std::vector<double> arcLength(const std::vector<Point2D>& pts);

// Extract coordinates
std::vector<double> getX(const std::vector<Point2D>& pts);
std::vector<double> getY(const std::vector<Point2D>& pts);

// 1D cubic spline
class CubicSpline1D {
public:
    void build(const std::vector<double>& s, const std::vector<double>& f);
    double eval(double sQuery) const;

private:
    std::size_t findInterval(double sQuery) const;
    double evalInterval(std::size_t i, double sQuery) const;

    std::vector<double> s_, a_, b_, c_, d_;
};

// q = 1,2,3 Lagrange edge-node locations
std::vector<double> lagrangePoints(int q);

// 2D spline wrapper
struct ParametricSpline2D {
    CubicSpline1D sx, sy;
    Point2D eval(double s) const;
};

// Recover spline parameter of a boundary point
double findS(const Point2D& p, const ParametricSpline2D& spline, double s0, double s1);

// Construct curved edge nodes
std::vector<Point2D> curvedEdgeNodes(
    const Point2D& A,
    const Point2D& B,
    const ParametricSpline2D& spline,
    double s0,
    double s1,
    int q
);

// Read nodes, curve-node pairs from .gri
GriMesh readGri(const std::string& filename);

// Build curved edges for one named boundary
std::vector<CurvedBoundaryEdge> buildCurvedEdges(
    const std::string& curveName,
    const std::vector<BoundaryPair>& edges,
    const std::vector<Point2D>& meshNodes,
    const ParametricSpline2D& spline,
    double s0,
    double s1,
    int q
);

// Write curved edges to file
void writeCurvedEdges(
    const std::string& filename,
    const std::vector<CurvedBoundaryEdge>& edges
);

#endif