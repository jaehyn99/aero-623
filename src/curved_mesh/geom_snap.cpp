#include "curved_mesh/geom_snap.h"
// #include "geom_snap.h"
#include <fstream>
#include <iostream>
#include <cmath>

// Read blade geometry points
std::vector<Point2D> readBladePoints(const std::string& filename) {
    std::ifstream in(filename);
    std::vector<Point2D> pts;
    double x, y;
    while (in >> x >> y) {
        pts.push_back({x, y});
    }
    return pts;
}

// Offset coordinates
std::vector<Point2D> offset(const std::vector<Point2D>& pts, double dx, double dy) {
    std::vector<Point2D> shifted = pts;
    for (auto& p : shifted) {
        p.x += dx;
        p.y += dy;
    }
    return shifted;
}

// Arc-length parameterization
std::vector<double> arcLength(const std::vector<Point2D>& pts) {
    std::vector<double> s(pts.size(), 0.0);
    for (std::size_t i = 1; i < pts.size(); ++i) {
        double dx = pts[i].x - pts[i-1].x;
        double dy = pts[i].y - pts[i-1].y;
        s[i] = s[i-1] + std::sqrt(dx*dx + dy*dy);
    }
    return s;
}

std::vector<double> getX(const std::vector<Point2D>& pts) {
    std::vector<double> x(pts.size());
    for (std::size_t i = 0; i < pts.size(); ++i) x[i] = pts[i].x;
    return x;
}
std::vector<double> getY(const std::vector<Point2D>& pts) {
    std::vector<double> y(pts.size());
    for (std::size_t i = 0; i < pts.size(); ++i) y[i] = pts[i].y;
    return y;
}

void CubicSpline1D::build(const std::vector<double>& s, const std::vector<double>& f) {
    const std::size_t n = s.size();
    s_ = s;
    a_ = f;

    std::vector<double> h(n - 1);
    for (std::size_t i = 0; i < n - 1; ++i) {
        h[i] = s_[i + 1] - s_[i];
    }

    std::vector<double> alpha(n, 0.0);
    for (std::size_t i = 1; i < n - 1; ++i) {
        alpha[i] = (3.0 / h[i]) * (a_[i + 1] - a_[i])
                 - (3.0 / h[i - 1]) * (a_[i] - a_[i - 1]);
    }

    std::vector<double> l(n), mu(n), z(n);
    c_.assign(n, 0.0);
    b_.assign(n - 1, 0.0);
    d_.assign(n - 1, 0.0);

    l[0] = 1.0;
    mu[0] = 0.0;
    z[0] = 0.0;

    for (std::size_t i = 1; i < n - 1; ++i) {
        l[i] = 2.0 * (s_[i + 1] - s_[i - 1]) - h[i - 1] * mu[i - 1];
        mu[i] = h[i] / l[i];
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
    }

    l[n - 1] = 1.0;
    z[n - 1] = 0.0;
    c_[n - 1] = 0.0;

    for (int j = static_cast<int>(n) - 2; j >= 0; --j) {
        c_[j] = z[j] - mu[j] * c_[j + 1];
        b_[j] = (a_[j + 1] - a_[j]) / h[j]
              - h[j] * (2.0 * c_[j] + c_[j + 1]) / 3.0;
        d_[j] = (c_[j + 1] - c_[j]) / (3.0 * h[j]);
    }
}

double CubicSpline1D::eval(double sQuery) const {
    if (sQuery <= s_.front()) return evalInterval(0, sQuery);
    if (sQuery >= s_.back())  return evalInterval(s_.size() - 2, sQuery);

    std::size_t i = findInterval(sQuery);
    return evalInterval(i, sQuery);
}

std::size_t CubicSpline1D::findInterval(double sQuery) const {
    for (std::size_t i = 0; i < s_.size() - 1; ++i) {
        if (sQuery >= s_[i] && sQuery <= s_[i + 1]) return i;
    }
    return s_.size() - 2;
}

double CubicSpline1D::evalInterval(std::size_t i, double sQuery) const {
    double ds = sQuery - s_[i];
    return a_[i] + b_[i]*ds + c_[i]*ds*ds + d_[i]*ds*ds*ds;
}

Point2D ParametricSpline2D::eval(double s) const {
    return {sx.eval(s), sy.eval(s)};
}

// q = 1,2,3 Lagrange edge-node locations
std::vector<double> lagrangePoints(int q) {
    if (q == 1) return {0.0, 1.0};
    if (q == 2) return {0.0, 0.5, 1.0};
    return {0.0, 1.0/3.0, 2.0/3.0, 1.0};
}

// recover spline parameter of a boundary point
double findS(const Point2D& p, const ParametricSpline2D& spline, double s0, double s1) {
    double sBest = s0;
    double d2Best = 1e300;

    for (int i = 0; i <= 4000; ++i) {
        double t = static_cast<double>(i) / 4000.0;
        double s = s0 + t * (s1 - s0);

        Point2D q = spline.eval(s);
        double dx = q.x - p.x;
        double dy = q.y - p.y;
        double d2 = dx*dx + dy*dy;

        if (d2 < d2Best) {
            d2Best = d2;
            sBest = s;
        }
    }
    return sBest;
}

std::vector<Point2D> curvedEdgeNodes(
    const Point2D& A,
    const Point2D& B,
    const ParametricSpline2D& spline,
    double s0,
    double s1,
    int q)
{
    double sA = findS(A, spline, s0, s1);
    double sB = findS(B, spline, s0, s1);

    auto xi = lagrangePoints(q);
    std::vector<Point2D> nodes;
    nodes.reserve(xi.size());

    for (double r : xi) {
        double s = (1.0 - r) * sA + r * sB;
        nodes.push_back(spline.eval(s));
    }
    return nodes;
}

// Read nodes, curve-node pairs from .gri
GriMesh readGri(const std::string& filename) {
    std::ifstream in(filename);
    GriMesh mesh;

    int Nn, Ne, dim;
    in >> Nn >> Ne >> dim;

    mesh.nodes.resize(Nn);
    for (int i = 0; i < Nn; ++i) {
        in >> mesh.nodes[i].x >> mesh.nodes[i].y;
    }

    int nCurves;
    in >> nCurves;

    for (int c = 0; c < nCurves; ++c) {
        int nEdges, nodesPerEdge;
        std::string name;
        in >> nEdges >> nodesPerEdge >> name;

        std::vector<BoundaryPair> pairs(nEdges);
        for (int i = 0; i < nEdges; ++i) {
            in >> pairs[i].n1 >> pairs[i].n2;
        }
        mesh.curves[name] = pairs;
    }

    return mesh;
}

std::vector<CurvedBoundaryEdge> buildCurvedEdges(
    const std::string& curveName,
    const std::vector<BoundaryPair>& edges,
    const std::vector<Point2D>& meshNodes,
    const ParametricSpline2D& spline,
    double s0,
    double s1,
    int q)
{
    std::vector<CurvedBoundaryEdge> curved;
    curved.reserve(edges.size());

    for (const auto& e : edges) {
        Point2D A = meshNodes[e.n1 - 1];
        Point2D B = meshNodes[e.n2 - 1];

        CurvedBoundaryEdge ce;
        ce.curveName = curveName;
        ce.n1 = e.n1;
        ce.n2 = e.n2;
        ce.qnodes = curvedEdgeNodes(A, B, spline, s0, s1, q);

        curved.push_back(ce);
    }
    return curved;
}


void writeCurvedEdges(const std::string& filename,
                      const std::vector<CurvedBoundaryEdge>& edges)
{
    std::ofstream out(filename);

    for (const auto& e : edges) {
        out << e.curveName << " " << e.n1 << " " << e.n2 << "\n";
        for (const auto& p : e.qnodes) {
            out << p.x << " " << p.y << "\n";
        }
        out << "\n";
    }
}