#include "FirstorderEuler.h"
#include "mesh/TriangularMesh.h"
#include "solver/BoundaryConditions.h"
#include "solver/boundaryFlux.hpp"
#include "solver/hlleFluxFO.hpp"
#include "solver/inletFlux.hpp"
#include "solver/numericalFlux.hpp"
#include "solver/outletFlux.hpp"
#include "solver/roeFlux.hpp"
#include "solver/wallFlux.hpp"

#include <algorithm>
#include <cmath>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <unordered_map>
#include <string>
#include <sstream>
#include <stdexcept>

namespace {

using Vec2 = FirstorderEuler::Vec2;
using BoundaryKind = FirstorderEuler::BoundaryKind;

inline double dot2(const Vec2& a, const Vec2& b) { return a[0]*b[0] + a[1]*b[1]; }
inline Vec2 neg2(const Vec2& a) { return Vec2{-a[0], -a[1]}; }
}

Vec2 FirstorderEuler::cellCentroid(std::size_t ei) const {
    const auto& e = elements_.at(ei);
    const auto& a = nodes_.at(e[0]);
    const auto& b = nodes_.at(e[1]);
    const auto& c = nodes_.at(e[2]);
    return Vec2{ (a[0]+b[0]+c[0])/3.0, (a[1]+b[1]+c[1])/3.0 };
}

namespace {

std::string toLower(std::string s) {
    for (char& c : s) {
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    }
    return s;
}

int parseCurveId(const std::string& title) {
    std::string digits;
    for (char c : title) {
        if (std::isdigit(static_cast<unsigned char>(c))) {
            digits.push_back(c);
        }
    }
    if (digits.empty()) {
        return -1;
    }
    return std::stoi(digits);
}

BoundaryKind boundaryKindFromTitle(const std::string& title,
                                   const FirstorderEuler::SolverConfig& config) {
    const int curveID = parseCurveId(title);

    if (curveID == 2 || curveID == 4 || curveID == 6 || curveID == 8) {
        return BoundaryKind::Periodic;
    }
    if (curveID == 1 || curveID == 5) {
        return BoundaryKind::WallSlip;
    }
    if (curveID == config.outflowCurve) {
        return BoundaryKind::OutflowSubsonic;
    }
    if (curveID == config.inflowCurve) {
        return BoundaryKind::InflowSteady;
    }

    if (curveID == 3 || curveID == 7) {
        return (curveID == 3) ? BoundaryKind::InflowSteady : BoundaryKind::OutflowSubsonic;
    }

    return BoundaryKind::InflowSteady;
}

EulerBoundaryConditions makeBoundaryConditions(const FirstorderEuler::SolverConfig& config) {
    EulerBoundaryConditions::Config bc;
    bc.gamma = config.gamma;
    bc.gasConstant = config.gasConstant;
    bc.alpha = config.alpha;
    bc.pout = config.pout;
    bc.rho0 = config.rho0;
    bc.inflowCurve = config.inflowCurve;
    bc.outflowCurve = config.outflowCurve;

    bc.Tt = (config.a0 * config.a0) / (config.gamma * config.gasConstant);
    bc.pt = config.rho0 * config.gasConstant * bc.Tt;
    bc.Vrot = config.a0;
    return EulerBoundaryConditions(bc);
}

EulerBoundaryConditions::Type toBoundaryType(BoundaryKind kind) {
    switch (kind) {
    case BoundaryKind::InflowSteady:
        return EulerBoundaryConditions::Type::InflowSteady;
    case BoundaryKind::InflowUnsteady:
        return EulerBoundaryConditions::Type::InflowUnsteady;
    case BoundaryKind::OutflowSubsonic:
        return EulerBoundaryConditions::Type::OutflowSubsonic;
    case BoundaryKind::WallSlip:
        return EulerBoundaryConditions::Type::WallSlip;
    case BoundaryKind::Periodic:
        return EulerBoundaryConditions::Type::Periodic;
    }
    return EulerBoundaryConditions::Type::InflowSteady;
}


FirstorderEuler::Vec2 normalized(const FirstorderEuler::Vec2& n) {
    const double mag = std::sqrt(n[0] * n[0] + n[1] * n[1]);
    if (mag <= 1e-14) {
        throw std::runtime_error("Encountered near-zero face normal.");
    }
    return {n[0] / mag, n[1] / mag};
}

double spectralRadius(const FirstorderEuler::Conserved& U,
                      const FirstorderEuler::Vec2& n,
                      double gamma);

Eigen::Vector4d toEigen(const FirstorderEuler::Conserved& u) {
    return Eigen::Vector4d(u[0], u[1], u[2], u[3]);
}

FirstorderEuler::Conserved fromEigen(const Eigen::Vector4d& u) {
    return {u(0), u(1), u(2), u(3)};
}

bool hasReasonablePressureAndInternalEnergy(const FirstorderEuler::Conserved& Unew,
                                            const FirstorderEuler::Conserved& Uold,
                                            double gamma) {
    constexpr double rhoFloor = 1e-10;
    constexpr double pFloor = 1e-10;
    constexpr double minPressureRatio = 0.7;

    const double rhoNew = Unew[0];
    if (!std::isfinite(rhoNew) || rhoNew <= rhoFloor) {
        return false;
    }

    const auto pressureOf = [&](const FirstorderEuler::Conserved& U) {
        const double rho = std::max(1e-14, U[0]);
        const double u = U[1] / rho;
        const double v = U[2] / rho;
        return (gamma - 1.0) * (U[3] - 0.5 * rho * (u * u + v * v));
    };

    const double pOld = pressureOf(Uold);
    const double pNew = pressureOf(Unew);
    if (!std::isfinite(pOld) || !std::isfinite(pNew) || pNew <= pFloor) {
        return false;
    }

    const double minAllowed = minPressureRatio * std::max(pFloor, pOld);
    if (pNew < minAllowed) {
        return false;
    }

    const double rho = rhoNew;
    const double u = Unew[1] / rho;
    const double v = Unew[2] / rho;
    const double eInt = Unew[3] - 0.5 * rho * (u * u + v * v);
    return std::isfinite(eInt) && eInt > pFloor / (gamma - 1.0);
}

void enforcePhysicalState(FirstorderEuler::Conserved& U, double gamma) {
    constexpr double rhoFloor = 1e-10;
    constexpr double pFloor = 1e-10;

    double rho = std::isfinite(U[0]) ? U[0] : rhoFloor;
    rho = std::max(rhoFloor, rho);

    double rhou = std::isfinite(U[1]) ? U[1] : 0.0;
    double rhov = std::isfinite(U[2]) ? U[2] : 0.0;
    double rhoE = std::isfinite(U[3]) ? U[3] : pFloor / (gamma - 1.0);

    const double u = rhou / rho;
    const double v = rhov / rho;
    const double kinetic = 0.5 * rho * (u * u + v * v);
    const double p = (gamma - 1.0) * (rhoE - kinetic);
    if (!std::isfinite(p) || p < pFloor) {
        rhoE = pFloor / (gamma - 1.0) + kinetic;
    }

    U[0] = rho;
    U[1] = rhou;
    U[2] = rhov;
    U[3] = rhoE;
}

std::unique_ptr<numericalFlux> makeFlux(const std::string& schemeName) {
    const std::string name = toLower(schemeName);
    if (name == "roe") {
        return std::make_unique<RoeFlux>();
    }
    if (name == "hlle") {
        return std::make_unique<HLLEFluxFO>();
    }
    throw std::runtime_error("Unknown flux scheme: " + schemeName + " (expected 'roe' or 'hlle').");
}

double spectralRadius(const FirstorderEuler::Conserved& U,
                      const FirstorderEuler::Vec2& n,
                      double gamma) {
    for (double x : U) {
        if (!std::isfinite(x)) {
            return 0.0;
        }
    }
    const double rho = std::max(1e-14, U[0]);
    const double u = U[1] / rho;
    const double v = U[2] / rho;
    const double p = std::max(1e-14, (gamma - 1.0) * (U[3] - 0.5 * rho * (u * u + v * v)));
    const double c = std::sqrt(gamma * p / rho);
    const double un = u * n[0] + v * n[1];
    const double s = std::abs(un) + c;
    return std::isfinite(s) ? s : 0.0;
}

} // namespace


FirstorderEuler::FirstorderEuler(MeshInputs inputs, SolverConfig config)
    : inputs_(std::move(inputs)), config_(std::move(config)) {}

void FirstorderEuler::loadInputs() {
    readMeshAndConnectivity();
    if (config_.validateMeshOnLoad) {
        validateLoadedArrays();
    }

    U_.assign(elements_.size(), Conserved{0.0, 0.0, 0.0, 0.0});
    residual_.assign(elements_.size(), Conserved{0.0, 0.0, 0.0, 0.0});

    if (config_.enableDebugPrints) {
        printBoundaryConditionSummary();
    }
}

void FirstorderEuler::readMeshAndConnectivity() {
    // Simpler flow: parse .gri once with TriangularMesh and derive element/face arrays from it.
    triMesh_ = std::make_unique<TriangularMesh>(inputs_.meshFile);

    nodes_.clear();
    nodes_.reserve(static_cast<std::size_t>(triMesh_->numNodes()));
    for (std::size_t i = 0; i < static_cast<std::size_t>(triMesh_->numNodes()); ++i) {
        const auto& p = triMesh_->node(i);
        nodes_.push_back({p[0], p[1]});
    }

    elements_.clear();
    elements_.reserve(static_cast<std::size_t>(triMesh_->numElems()));
    for (std::size_t i = 0; i < static_cast<std::size_t>(triMesh_->numElems()); ++i) {
        const auto& e = triMesh_->elem(i);
        elements_.push_back({static_cast<std::size_t>(e._pointID[0]),
                             static_cast<std::size_t>(e._pointID[1]),
                             static_cast<std::size_t>(e._pointID[2])});
    }

    area_.clear();
    area_.reserve(static_cast<std::size_t>(triMesh_->numElems()));
    for (std::size_t i = 0; i < static_cast<std::size_t>(triMesh_->numElems()); ++i) {
        const double ai = triMesh_->area(i);
        if (ai <= 0.0) {
            throw std::runtime_error("TriangularMesh area is non-positive.");
        }
        area_.push_back(ai);
    }

    buildFacesFromTriangularMesh();
    computePerimeterFromMesh();
}

void FirstorderEuler::buildFacesFromTriangularMesh() {
    if (!triMesh_) {
        throw std::runtime_error("TriangularMesh must be initialized before building face data.");
    }

    interiorFaces_.clear();
    boundaryFaces_.clear();
    periodicEdges_.clear();

    // Build face->owner-element mapping directly from element connectivity.
    // This avoids depending on TriangularMesh face owner IDs, so solver-side
    // assembly stays robust even if another branch has older mesh code.
    std::vector<std::array<int, 2>> faceOwners(static_cast<std::size_t>(triMesh_->numFaces()), {-1, -1});
    for (std::size_t elemID = 0; elemID < static_cast<std::size_t>(triMesh_->numElems()); ++elemID) {
        const auto& elem = triMesh_->elem(elemID);
        for (int fid : elem._faceID) {
            if (fid < 0 || static_cast<std::size_t>(fid) >= faceOwners.size()) {
                throw std::runtime_error("Element contains out-of-range face ID.");
            }
            auto& owners = faceOwners[static_cast<std::size_t>(fid)];
            if (owners[0] == -1) {
                owners[0] = static_cast<int>(elemID);
            } else if (owners[1] == -1) {
                owners[1] = static_cast<int>(elemID);
            } else {
                throw std::runtime_error("Face has more than two owner elements.");
            }
        }
    }

    std::vector<int> periodicPartnerFace(static_cast<std::size_t>(triMesh_->numFaces()), -1);

    std::unordered_map<std::string, std::vector<std::size_t>> boundaryFacesByTitle;
    for (std::size_t fid = 0; fid < static_cast<std::size_t>(triMesh_->numFaces()); ++fid) {
        const auto& face = triMesh_->face(fid);
        if (face.isBoundaryFace()) {
            boundaryFacesByTitle[face._title].push_back(fid);
        }
    }

    const auto pairPeriodicTitles = [&](const std::string& titleA, const std::string& titleB) {
        const auto itA = boundaryFacesByTitle.find(titleA);
        const auto itB = boundaryFacesByTitle.find(titleB);
        if (itA == boundaryFacesByTitle.end() || itB == boundaryFacesByTitle.end()) {
            return;
        }
        auto facesA = itA->second;
        auto facesB = itB->second;
        if (facesA.size() != facesB.size()) {
            std::ostringstream oss;
            oss << "Periodic boundary mismatch: " << titleA << " has " << facesA.size()
                << " faces while " << titleB << " has " << facesB.size() << ".";
            throw std::runtime_error(oss.str());
        }

        const auto centerX = [&](std::size_t fid) {
            const auto& f = triMesh_->face(fid);
            const auto& p0 = triMesh_->node(static_cast<std::size_t>(f._pointID[0]));
            const auto& p1 = triMesh_->node(static_cast<std::size_t>(f._pointID[1]));
            return 0.5 * (p0[0] + p1[0]);
        };
        std::sort(facesA.begin(), facesA.end(), [&](std::size_t a, std::size_t b) { return centerX(a) < centerX(b); });
        std::sort(facesB.begin(), facesB.end(), [&](std::size_t a, std::size_t b) { return centerX(a) < centerX(b); });

        for (std::size_t i = 0; i < facesA.size(); ++i) {
            const std::size_t a = facesA[i];
            const std::size_t b = facesB[i];
            periodicPartnerFace[a] = static_cast<int>(b);
            periodicPartnerFace[b] = static_cast<int>(a);
        }
    };
    pairPeriodicTitles("Curve2", "Curve4");
    pairPeriodicTitles("Curve6", "Curve8");

    std::unordered_map<std::string, std::size_t> groupFromTitle;
    const auto findLocalFace = [&](std::size_t elemID,
                                   int faceID,
                                   std::size_t globalFaceID,
                                   const char* side) -> std::size_t {
        const auto& faceIDs = triMesh_->elem(elemID)._faceID;
        const auto it = std::find(faceIDs.begin(), faceIDs.end(), faceID);
        if (it == faceIDs.end()) {
            std::ostringstream oss;
            oss << "Face-to-element mapping error on " << side
                << " side: globalFaceID=" << globalFaceID
                << ", requestedFaceID=" << faceID
                << ", elemID=" << elemID
                << ", elemFaceIDs=[" << faceIDs[0] << "," << faceIDs[1] << "," << faceIDs[2] << "]";
            throw std::runtime_error(oss.str());
        }
        return static_cast<std::size_t>(it - faceIDs.begin());
    };



    for (std::size_t faceID = 0; faceID < static_cast<std::size_t>(triMesh_->numFaces()); ++faceID) {
        const auto& face = triMesh_->face(faceID);
        const auto& owners = faceOwners[faceID];

        const bool periodic = (periodicPartnerFace[faceID] != -1);
        // IMPORTANT: don't process both faces of the periodic pair
        if (periodic) {
            const int partner = periodicPartnerFace[faceID];
            if (partner < 0) throw std::runtime_error("Periodic face has no partner.");
            if (static_cast<int>(faceID) > partner) continue; // skip duplicate half
        }

        const bool boundary = face.isBoundaryFace() && !periodic;

        if (owners[0] < 0) {
            throw std::runtime_error("Face has invalid left element index.");
        }

        const std::size_t elemL = static_cast<std::size_t>(owners[0]);
        const std::size_t localFaceL = findLocalFace(elemL, static_cast<int>(faceID), faceID, "left");

        if (boundary) {
            const auto [it, inserted] = groupFromTitle.try_emplace(face._title, groupFromTitle.size());
            (void)inserted;

            BoundaryFace bf;
            bf.elem = elemL;
            bf.localFace = localFaceL;
            bf.boundaryGroup = it->second;
            bf.boundaryTitle = face._title;
            bf.normal = {triMesh_->normal(elemL, localFaceL)[0], triMesh_->normal(elemL, localFaceL)[1]};
            const auto& p0 = triMesh_->node(static_cast<std::size_t>(face._pointID[0]));
            const auto& p1 = triMesh_->node(static_cast<std::size_t>(face._pointID[1]));
            bf.center = {(p0[0] + p1[0]) * 0.5, (p0[1] + p1[1]) * 0.5};
            bf.length = face._length;
            boundaryFaces_.push_back(bf);
            continue;
        }

        std::size_t elemR = 0;
        std::size_t localFaceR = 0;
        if (periodic) {
            const std::size_t periodicFaceID = static_cast<std::size_t>(periodicPartnerFace[faceID]);
            if (periodicFaceID >= faceOwners.size() || faceOwners[periodicFaceID][0] < 0) {
                throw std::runtime_error("Periodic face partner has invalid owner element.");
            }
            elemR = static_cast<std::size_t>(faceOwners[periodicFaceID][0]);
            localFaceR = findLocalFace(elemR, static_cast<int>(periodicFaceID), periodicFaceID, "right-periodic");
        } else {
            if (owners[1] < 0) {
                throw std::runtime_error("Interior face has invalid right element index.");
            }
            elemR = static_cast<std::size_t>(owners[1]);
            localFaceR = findLocalFace(elemR, static_cast<int>(faceID), faceID, "right-interior");
        }

        InteriorFace inf;
        inf.elemL = elemL;
        inf.faceL = localFaceL;
        inf.elemR = elemR;
        inf.faceR = localFaceR;
        inf.normal = {triMesh_->normal(elemL, localFaceL)[0], triMesh_->normal(elemL, localFaceL)[1]};
        const auto& p0 = triMesh_->node(static_cast<std::size_t>(face._pointID[0]));
        const auto& p1 = triMesh_->node(static_cast<std::size_t>(face._pointID[1]));
        inf.center = {(p0[0] + p1[0]) * 0.5, (p0[1] + p1[1]) * 0.5};
        inf.length = face._length;
        interiorFaces_.push_back(inf);

        if (periodic) {
            periodicEdges_.push_back({elemL, elemR});
        }

        //Debug
        if (config_.enableDebugPrints) {
            std::cout << "[debug] built faces: interior=" << interiorFaces_.size()
              << " boundary=" << boundaryFaces_.size()
              << " periodicPairs=" << periodicEdges_.size() << "\n";
        }
    }
}

void FirstorderEuler::computePerimeterFromMesh() {
    if (!triMesh_) {
        throw std::runtime_error("TriangularMesh must be initialized before computing perimeter.");
    }

    perimeter_.assign(elements_.size(), 0.0);
    for (std::size_t elemID = 0; elemID < static_cast<std::size_t>(triMesh_->numElems()); ++elemID) {
        const auto& elem = triMesh_->elem(elemID);
        perimeter_[elemID] = triMesh_->length(static_cast<std::size_t>(elem._faceID[0]))
            + triMesh_->length(static_cast<std::size_t>(elem._faceID[1]))
            + triMesh_->length(static_cast<std::size_t>(elem._faceID[2]));
    }
}

void FirstorderEuler::validateLoadedArrays() const {
    if (nodes_.empty() || elements_.empty()) {
        throw std::runtime_error("Mesh must contain nodes and elements.");
    }
    if (area_.size() != elements_.size()) {
        throw std::runtime_error("TriangularMesh area count must equal number of elements.");
    }

    for (const auto& f : interiorFaces_) {
        if (f.elemL >= elements_.size() || f.elemR >= elements_.size()) {
            throw std::runtime_error("I2E element index out of range.");
        }
        if (f.faceL > 2 || f.faceR > 2) {
            throw std::runtime_error("I2E local face index must be in {1,2,3} (0-based {0,1,2}).");
        }
    }

    for (const auto& f : boundaryFaces_) {
        if (f.elem >= elements_.size()) {
            throw std::runtime_error("B2E element index out of range.");
        }
        if (f.localFace > 2) {
            throw std::runtime_error("B2E local face index must be in {1,2,3} (0-based {0,1,2}).");
        }
    }
}

void FirstorderEuler::initUniformState() {
    const double M = config_.initialMach;
    const double g = config_.gamma;
    const double R = config_.gasConstant;
    const double p0 = config_.rho0 * config_.a0 * config_.a0 / g;

    const double T0 = (config_.a0 * config_.a0) / (g * R);
    const double T = T0 / (1.0 + 0.5 * (g - 1.0) * M * M);
    const double p = p0 * std::pow(T / T0, g / (g - 1.0));

    const double u = M * config_.a0 * std::cos(config_.alpha);
    const double v = M * config_.a0 * std::sin(config_.alpha);

    const double rho = p / (R * T);
    const double rhoE = p / (g - 1.0) + 0.5 * rho * (u * u + v * v);

    for (auto& Ui : U_) {
        Ui[0] = rho;
        Ui[1] = rho * u;
        Ui[2] = rho * v;
        Ui[3] = rhoE;
    }
}

void FirstorderEuler::advanceToConvergedOrFinalTime() {
    advance(true);
}

void FirstorderEuler::runSteadyGlobal() {
    config_.localTimeStepping = false;
    resetMarchState();
    advance(false);
}

void FirstorderEuler::runSteadyLocal() {
    config_.localTimeStepping = true;
    resetMarchState();
    advance(false);
}

void FirstorderEuler::runUnsteadyGlobal() {
    config_.localTimeStepping = false;
    resetMarchState();
    advance(true);
}

void FirstorderEuler::resetMarchState() {
    time_ = 0.0;
    iteration_ = 0;
}

void FirstorderEuler::writeSolutionVtk(const std::string& filePath) const {
    if (nodes_.empty() || elements_.empty() || U_.size() != elements_.size()) {
        throw std::runtime_error("Cannot export VTK: mesh/state arrays are not initialized.");
    }

    const double g = config_.gamma;
    const double MInf = std::max(1e-8, config_.initialMach);
    const double pInf = (config_.rho0 * config_.a0 * config_.a0 / g)
                      / std::pow(1.0 + 0.5 * (g - 1.0) * MInf * MInf, g / (g - 1.0));
    const double qInf = std::max(1e-12, 0.5 * pInf * g * MInf * MInf);

    std::ofstream out(filePath);
    if (!out) {
        throw std::runtime_error("Failed to open VTK output file: " + filePath);
    }

    out << "# vtk DataFile Version 3.0\n";
    out << "FirstorderEuler solution\n";
    out << "ASCII\n";
    out << "DATASET UNSTRUCTURED_GRID\n";

    out << "POINTS " << nodes_.size() << " float\n";
    out << std::setprecision(10);
    for (const auto& n : nodes_) {
        out << n[0] << " " << n[1] << " 0.0\n";
    }

    out << "CELLS " << elements_.size() << " " << (elements_.size() * 4) << "\n";
    for (const auto& e : elements_) {
        out << "3 " << e[0] << " " << e[1] << " " << e[2] << "\n";
    }

    out << "CELL_TYPES " << elements_.size() << "\n";
    for (std::size_t i = 0; i < elements_.size(); ++i) {
        out << "5\n";
    }

    out << "CELL_DATA " << elements_.size() << "\n";

    out << "SCALARS pressure float 1\n";
    out << "LOOKUP_TABLE default\n";
    for (const auto& U : U_) {
        const double rho = std::max(1e-14, U[0]);
        const double u = U[1] / rho;
        const double v = U[2] / rho;
        const double p = std::max(1e-14, (g - 1.0) * (U[3] - 0.5 * rho * (u * u + v * v)));
        out << p << "\n";
    }

    out << "SCALARS mach float 1\n";
    out << "LOOKUP_TABLE default\n";
    for (const auto& U : U_) {
        const double rho = std::max(1e-14, U[0]);
        const double u = U[1] / rho;
        const double v = U[2] / rho;
        const double p = std::max(1e-14, (g - 1.0) * (U[3] - 0.5 * rho * (u * u + v * v)));
        const double a = std::sqrt(g * p / rho);
        const double V = std::sqrt(u * u + v * v);
        out << (V / std::max(1e-14, a)) << "\n";
    }

    out << "SCALARS cp float 1\n";
    out << "LOOKUP_TABLE default\n";
    for (const auto& U : U_) {
        const double rho = std::max(1e-14, U[0]);
        const double u = U[1] / rho;
        const double v = U[2] / rho;
        const double p = std::max(1e-14, (g - 1.0) * (U[3] - 0.5 * rho * (u * u + v * v)));
        out << ((p - pInf) / qInf) << "\n";
    }
}

void FirstorderEuler::advance(bool stopByTime) {
    std::cout << "    it           t          dt        ||R||2\n";

    while (iteration_ < config_.maxIterations && (!stopByTime || time_ < config_.finalTime)) {
        auto edges = computeEdgeFluxesAndWaveSpeeds();
        assembleResidualFromEdgeFluxes(edges);

        double dtUsed = 0.0;
        std::vector<double> dtLocalBuffer;
        const std::vector<double>* dtLocalPtr = nullptr;
        if (config_.localTimeStepping) {
            dtLocalBuffer = computeLocalDtFromWaveSpeeds(edges);
            dtLocalPtr = &dtLocalBuffer;
            dtUsed = dtLocalBuffer.empty() ? 0.0 : *std::min_element(dtLocalBuffer.begin(), dtLocalBuffer.end());
            updateStateLocalDt(dtLocalBuffer);
        } else {
            dtUsed = computeGlobalDtFromWaveSpeeds(edges);
            updateStateGlobalDt(dtUsed);
        }

        if (!std::isfinite(dtUsed) || dtUsed <= 0.0) {
            throw std::runtime_error("Computed non-finite/invalid dt.");
        }

        const double normR = l2Norm(residual_);
        if (!std::isfinite(normR)) {
            const std::string dumpName = config_.outputPrefix + "_nan_iter" + std::to_string(iteration_) + ".vtk";
            writeSolutionVtk(dumpName);
            throw std::runtime_error("Residual norm became non-finite (NaN/Inf).\n"
                                     "Wrote debug VTK: " + dumpName + "\n"
                                     "Try smaller CFL or verify BC/mesh consistency.");
        }
        if (iteration_ % 10 == 0 || (stopByTime && time_ >= config_.finalTime)) {
            std::cout << iteration_ << " " << time_ << " " << dtUsed << " " << normR << "\n";
        }

        if (config_.enableDebugPrints && (iteration_ % std::max<std::size_t>(1, config_.debugEvery) == 0)) {
            printIterationDiagnostics(edges, dtLocalPtr, dtUsed, normR);
        }

        ++iteration_;
        if (stopByTime) {
            time_ += dtUsed;
        }

        if (config_.saveEvery > 0 && iteration_ % config_.saveEvery == 0) {
            const std::string dumpName = config_.outputPrefix + "_iter" + std::to_string(iteration_) + ".vtk";
            writeSolutionVtk(dumpName);
        }

        if (normR < config_.residualTolerance) {
            break;
        }
    }
}

std::vector<FirstorderEuler::EdgeFluxContribution> FirstorderEuler::computeEdgeFluxesAndWaveSpeeds() const {
    std::vector<EdgeFluxContribution> edges;
    edges.reserve(interiorFaces_.size() + boundaryFaces_.size());

    const auto flux = makeFlux(config_.fluxScheme);

    for (const auto& f : interiorFaces_) {
        const Conserved& UL = U_.at(f.elemL);
        const Conserved& UR = U_.at(f.elemR);

        // --- (1) enforce n points from elemL -> elemR ---
        Vec2 nFace = f.normal;
        const Vec2 cL = cellCentroid(f.elemL);
        const Vec2 cR = cellCentroid(f.elemR);
        const Vec2 dLR{ cR[0] - cL[0], cR[1] - cL[1] };
        if (dot2(nFace, dLR) < 0.0) {
            nFace = neg2(nFace);
        }

        const Vec2 nUnit = normalized(nFace);
        const Eigen::Vector2d n(nUnit[0], nUnit[1]);

        const Conserved F = fromEigen((*flux)(toEigen(UL), toEigen(UR), config_.gamma, n));

        EdgeFluxContribution edge;
        edge.ownerElem = f.elemL;
        edge.neighborElem = f.elemR;
        edge.normal = nFace;          // store the *oriented* normal
        edge.edgeLength = f.length;
        edge.flux = F;
        edge.spectralRadius = std::max(spectralRadius(UL, nUnit, config_.gamma),
                                    spectralRadius(UR, nUnit, config_.gamma));
        edges.push_back(edge);
    }



    for (const auto& f : boundaryFaces_) {
        const Conserved& UL = U_.at(f.elem);
        const BoundaryKind kind = boundaryKindFromTitle(f.boundaryTitle, config_);
        const EulerBoundaryConditions bcModel = makeBoundaryConditions(config_);

        Vec2 nFace = f.normal;
        const Vec2 c = cellCentroid(f.elem);
        const Vec2 toFace{ f.center[0] - c[0], f.center[1] - c[1] };
        if (dot2(nFace, toFace) < 0.0) nFace = neg2(nFace);

        BoundaryFace fFixed = f;
        fFixed.normal = nFace;

        const Vec2 nUnit = normalized(nFace);
        const Eigen::Vector2d n(nUnit[0], nUnit[1]);

        EulerBoundaryConditions::Context ctx;
        ctx.time = time_;
        ctx.faceCenter = fFixed.center;
        const Conserved UR = bcModel.boundaryState(toBoundaryType(kind), UL, nUnit, ctx);

        const Conserved F = fromEigen((*flux)(toEigen(UL), toEigen(UR), config_.gamma, n));

        EdgeFluxContribution edge;
        edge.ownerElem = f.elem;
        edge.neighborElem = f.elem;
        edge.normal = nFace;
        edge.edgeLength = f.length;
        edge.flux = F;
        edge.spectralRadius = std::max(
            spectralRadius(UL, nUnit, config_.gamma),
            spectralRadius(UR, nUnit, config_.gamma)
        );

        edges.push_back(edge);
    }

    return edges;
}

FirstorderEuler::Conserved FirstorderEuler::computeBoundaryFluxFromModules(
    const BoundaryFace& f,
    const Conserved& UL,
    BoundaryKind kind) const {
    const Eigen::Vector4d up = toEigen(UL);
    const Vec2 nUnit = normalized(f.normal);
    const Eigen::Vector2d n(nUnit[0], nUnit[1]);

    std::unique_ptr<boundaryFlux> faceFlux;
    switch (kind) {
    case BoundaryKind::InflowSteady:
    case BoundaryKind::InflowUnsteady: {
        const bool transient = (kind == BoundaryKind::InflowUnsteady);
        faceFlux = std::make_unique<inletFlux>(config_.rho0, config_.a0, config_.alpha, time_, transient);
        break;
    }
    case BoundaryKind::OutflowSubsonic:
        faceFlux = std::make_unique<outletFlux>(config_.pout);
        break;
    case BoundaryKind::WallSlip:
        faceFlux = std::make_unique<wallFlux>();
        break;
    case BoundaryKind::Periodic:
        throw std::runtime_error("Periodic boundary was routed to boundary-flux evaluation.");
    }

    return fromEigen((*faceFlux)(up, config_.gamma, n));
}

void FirstorderEuler::assembleResidualFromEdgeFluxes(const std::vector<EdgeFluxContribution>& edges) {
    for (auto& R : residual_) {
        for (double& x : R) {
            x = 0.0;
        }
    }

    for (const auto& edge : edges) {
        auto& RL = residual_.at(edge.ownerElem);
        for (std::size_t k = 0; k < RL.size(); ++k) {
            RL[k] += edge.edgeLength * edge.flux[k];
        }

        if (edge.neighborElem != edge.ownerElem) {
            auto& RR = residual_.at(edge.neighborElem);
            for (std::size_t k = 0; k < RR.size(); ++k) {
                RR[k] -= edge.edgeLength * edge.flux[k];
            }
        }
    }
}

double FirstorderEuler::computeGlobalDtFromWaveSpeeds(const std::vector<EdgeFluxContribution>& edges) const {
    if (U_.empty()) {
        return 0.0;
    }

    if (edges.empty()) {
        return 1e-12;
    }

    // Consistent with PDF Eq. (3.4.3)-(3.4.5):
    // d_i = 2*A_i/P_i,
    // |s|_i = sum_e |s|_{i,e} * Delta l_{i,e} / P_i,
    // global dt = min_i (CFL * d_i / |s|_i).
    std::vector<double> speedLenSum(area_.size(), 0.0);
    for (const auto& edge : edges) {
        const double contrib = std::abs(edge.spectralRadius) * edge.edgeLength;
        if (!std::isfinite(contrib) || contrib <= 0.0) {
            continue;
        }
        speedLenSum[edge.ownerElem] += contrib;
        if (edge.neighborElem != edge.ownerElem) {
            speedLenSum[edge.neighborElem] += contrib;
        }
    }

    double dtMin = std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < area_.size(); ++i) {
        const double Pi = std::max(1e-12, perimeter_[i]);
        const double di = 2.0 * area_[i] / Pi;
        const double sbar = speedLenSum[i] / Pi;
        const double dti = config_.cfl * di / std::max(1e-12, sbar);
        if (std::isfinite(dti) && dti > 0.0) {
            dtMin = std::min(dtMin, dti);
        }
    }

    if (!std::isfinite(dtMin)) {
        return 1e-12;
    }
    return std::max(1e-12, dtMin);
}

std::vector<double> FirstorderEuler::computeLocalDtFromWaveSpeeds(const std::vector<EdgeFluxContribution>& edges) const {
    std::vector<double> dt(area_.size(), 1e-3);
    if (edges.empty()) {
        return dt;
    }

    std::vector<double> edgeSpeedWeighted(area_.size(), 0.0);
    for (const auto& edge : edges) {
        const double s = std::abs(edge.spectralRadius) * edge.edgeLength;
        if (!std::isfinite(s) || s <= 0.0) {
            continue;
        }
        edgeSpeedWeighted[edge.ownerElem] += s;
        if (edge.neighborElem != edge.ownerElem) {
            edgeSpeedWeighted[edge.neighborElem] += s;
        }
    }

    for (std::size_t i = 0; i < area_.size(); ++i) {
        const double dchar = 2.0 * area_[i] / std::max(1e-12, perimeter_[i]);
        const double dtMax = std::max(1e-8, config_.cfl * dchar);
        const double dtRaw = 2.0 * area_[i] * config_.cfl / std::max(1e-12, edgeSpeedWeighted[i]);
        dt[i] = std::max(1e-12, std::min(dtRaw, dtMax));
    }
    return dt;
}

void FirstorderEuler::updateStateGlobalDt(double dt) {
    const int maxCuts = 10;
    constexpr double maxRelativeEnergyChange = 0.25;

    for (std::size_t i = 0; i < U_.size(); ++i) {
        double dtUse = dt;
        Conserved Uold = U_[i];

        for (int cut = 0; cut <= maxCuts; ++cut) {
            Conserved Utry = Uold;
            const double scale = dtUse / area_[i];
            for (std::size_t k = 0; k < 4; ++k) {
                Utry[k] -= scale * residual_[i][k];
            }

            const double oldE = std::max(1e-12, std::abs(Uold[3]));
            const double relEnergyChange = std::abs(Utry[3] - Uold[3]) / oldE;
            const bool energyStepOK = std::isfinite(relEnergyChange) && relEnergyChange <= maxRelativeEnergyChange;
            const bool pressureStepOK = hasReasonablePressureAndInternalEnergy(Utry, Uold, config_.gamma);

            if (energyStepOK && pressureStepOK) {
                U_[i] = Utry;
                break;
            }

            if (cut == maxCuts) {
                U_[i] = Utry;
                enforcePhysicalState(U_[i], config_.gamma);
            } else {
                dtUse *= 0.5;
            }
        }
    }
}


void FirstorderEuler::updateStateLocalDt(const std::vector<double>& dtLocal) {
    const int maxCuts = 10;
    constexpr double maxRelativeEnergyChange = 0.25;

    for (std::size_t i = 0; i < U_.size(); ++i) {
        double dtUse = dtLocal[i];
        const Conserved Uold = U_[i];

        for (int cut = 0; cut <= maxCuts; ++cut) {
            Conserved Utry = Uold;
            const double scale = dtUse / area_[i];
            for (std::size_t k = 0; k < Utry.size(); ++k) {
                Utry[k] -= scale * residual_[i][k];
            }

            const double oldE = std::max(1e-12, std::abs(Uold[3]));
            const double relEnergyChange = std::abs(Utry[3] - Uold[3]) / oldE;
            const bool energyStepOK = std::isfinite(relEnergyChange) && relEnergyChange <= maxRelativeEnergyChange;
            const bool pressureStepOK = hasReasonablePressureAndInternalEnergy(Utry, Uold, config_.gamma);
            if (energyStepOK && pressureStepOK) {
                U_[i] = Utry;
                break;
            }

            if (cut == maxCuts) {
                U_[i] = Utry;
                enforcePhysicalState(U_[i], config_.gamma);
            } else {
                dtUse *= 0.5;
            }
        }
    }
}


double FirstorderEuler::cellPressure(const Conserved& U) const {
    const double rho = std::max(1e-14, U[0]);
    const double u = U[1] / rho;
    const double v = U[2] / rho;
    return (config_.gamma - 1.0) * (U[3] - 0.5 * rho * (u * u + v * v));
}

void FirstorderEuler::printBoundaryConditionSummary() const {
    std::size_t nInflow = 0, nOutflow = 0, nWall = 0, nPeriodic = 0;
    std::unordered_map<std::string, std::size_t> perTitle;

    for (const auto& f : boundaryFaces_) {
        perTitle[f.boundaryTitle] += 1;
        switch (boundaryKindFromTitle(f.boundaryTitle, config_)) {
        case BoundaryKind::InflowSteady:
        case BoundaryKind::InflowUnsteady:
            ++nInflow;
            break;
        case BoundaryKind::OutflowSubsonic:
            ++nOutflow;
            break;
        case BoundaryKind::WallSlip:
            ++nWall;
            break;
        case BoundaryKind::Periodic:
            ++nPeriodic;
            break;
        }
    }

    std::cout << "[debug] BC summary: boundaryFaces=" << boundaryFaces_.size()
              << " inflow=" << nInflow
              << " outflow=" << nOutflow
              << " wall=" << nWall
              << " periodic=" << nPeriodic
              << " periodicPairs=" << periodicEdges_.size() << "\n";
    for (const auto& kv : perTitle) {
        std::cout << "[debug]   title='" << kv.first << "' faces=" << kv.second << "\n";
    }
}

void FirstorderEuler::printIterationDiagnostics(const std::vector<EdgeFluxContribution>& edges,
                                                const std::vector<double>* dtLocal,
                                                double dtUsed,
                                                double normR) const {
    const auto minmaxRho = std::minmax_element(U_.begin(), U_.end(),
        [](const Conserved& a, const Conserved& b) { return a[0] < b[0]; });

    double pMin = std::numeric_limits<double>::infinity();
    double pMax = -std::numeric_limits<double>::infinity();
    std::size_t pMinCell = 0, pMaxCell = 0;
    for (std::size_t i = 0; i < U_.size(); ++i) {
        const double p = cellPressure(U_[i]);
        if (p < pMin) { pMin = p; pMinCell = i; }
        if (p > pMax) { pMax = p; pMaxCell = i; }
    }

    double fluxMax = 0.0;
    std::size_t fluxMaxEdge = 0;
    for (std::size_t i = 0; i < edges.size(); ++i) {
        const double m = std::max({std::abs(edges[i].flux[0]), std::abs(edges[i].flux[1]),
                                   std::abs(edges[i].flux[2]), std::abs(edges[i].flux[3])});
        if (m > fluxMax) { fluxMax = m; fluxMaxEdge = i; }
    }

    Vec2 pMinCenter{0.0, 0.0};
    if (pMinCell < elements_.size()) {
        const auto& e = elements_[pMinCell];
        const auto& a = nodes_[e[0]];
        const auto& b = nodes_[e[1]];
        const auto& c = nodes_[e[2]];
        pMinCenter = {(a[0] + b[0] + c[0]) / 3.0, (a[1] + b[1] + c[1]) / 3.0};
    }

    Vec2 fluxMaxCenter{0.0, 0.0};
    if (fluxMaxEdge < edges.size()) {
        const auto& ef = edges[fluxMaxEdge];
        for (const auto& f : interiorFaces_) {
            if ((f.elemL == ef.ownerElem && f.elemR == ef.neighborElem)
                || (f.elemL == ef.neighborElem && f.elemR == ef.ownerElem)) {
                fluxMaxCenter = f.center;
                break;
            }
        }
    }

    double rMax = 0.0;
    std::size_t rMaxCell = 0, rMaxComp = 0;
    for (std::size_t i = 0; i < residual_.size(); ++i) {
        for (std::size_t k = 0; k < 4; ++k) {
            const double a = std::abs(residual_[i][k]);
            if (a > rMax) { rMax = a; rMaxCell = i; rMaxComp = k; }
        }
    }

    auto dumpCellFluxContrib = [&](std::size_t cellID) {
        std::cout << "[debug]   cell " << cellID << " face flux contributions:\n";
        for (std::size_t ei = 0; ei < edges.size(); ++ei) {
            const auto& ed = edges[ei];
            if (ed.ownerElem == cellID || ed.neighborElem == cellID) {
                // sign: owner adds +F, neighbor adds -F
                const double sgn = (ed.ownerElem == cellID) ? +1.0 : -1.0;
                const bool boundaryEdge = (ed.ownerElem == ed.neighborElem);

                const Conserved& Uowner = U_.at(ed.ownerElem);
                const double pOwner = cellPressure(Uowner);
                double pNeighbor = pOwner;
                if (!boundaryEdge) {
                    pNeighbor = cellPressure(U_.at(ed.neighborElem));
                }

                std::cout << "    edge#" << ei
                        << " sgn=" << sgn
                        << " owner=" << ed.ownerElem
                        << " neigh=" << ed.neighborElem
                        << " bnd=" << (boundaryEdge ? 1 : 0)
                        << " L=" << ed.edgeLength
                        << " spec=" << ed.spectralRadius
                        << " pOwner=" << pOwner
                        << " pNeigh=" << pNeighbor
                        << " F=[" << ed.flux[0] << "," << ed.flux[1] << "," << ed.flux[2] << "," << ed.flux[3] << "]\n";
            }
        }
    };
    dumpCellFluxContrib(pMinCell);
    dumpCellFluxContrib(pMaxCell);
    dumpCellFluxContrib(rMaxCell);

    // Targeted debug probe for TE instability tracking.
    // Focus pair requested in debugging session: cells 11 and 12, edge 57.
    const auto dumpCellPrimitiveAndScale = [&](std::size_t cellID) {
        if (cellID >= U_.size() || cellID >= area_.size() || cellID >= perimeter_.size()) {
            std::cout << "[debug]   probe cell " << cellID << " is out of range.\n";
            return;
        }
        const Conserved& Uc = U_[cellID];
        const double rho = std::max(1e-14, Uc[0]);
        const double u = Uc[1] / rho;
        const double v = Uc[2] / rho;
        const double p = cellPressure(Uc);
        const double rhoE = Uc[3];

        const double dchar = 2.0 * area_[cellID] / std::max(1e-12, perimeter_[cellID]);
        const double dtCell = (dtLocal != nullptr && cellID < dtLocal->size()) ? (*dtLocal)[cellID] : dtUsed;
        const double dtOverA = dtCell / std::max(1e-14, area_[cellID]);

        std::cout << "[debug]   probe cell=" << cellID
                  << " rho=" << rho
                  << " u=" << u
                  << " v=" << v
                  << " p=" << p
                  << " rhoE=" << rhoE
                  << " A=" << area_[cellID]
                  << " P=" << perimeter_[cellID]
                  << " dchar=" << dchar
                  << " dt=" << dtCell
                  << " dt/A=" << dtOverA
                  << " RE=" << residual_[cellID][3]
                  << "\n";
    };

    dumpCellPrimitiveAndScale(11);
    dumpCellPrimitiveAndScale(12);

    if (std::size_t(57) < edges.size()) {
        const auto& e57 = edges[57];
        const bool boundary57 = (e57.ownerElem == e57.neighborElem);

        double re57c11 = 0.0;
        double re57c12 = 0.0;
        if (e57.ownerElem == 11) {
            re57c11 += e57.edgeLength * e57.flux[3];
        }
        if (!boundary57 && e57.neighborElem == 11) {
            re57c11 -= e57.edgeLength * e57.flux[3];
        }
        if (e57.ownerElem == 12) {
            re57c12 += e57.edgeLength * e57.flux[3];
        }
        if (!boundary57 && e57.neighborElem == 12) {
            re57c12 -= e57.edgeLength * e57.flux[3];
        }

        std::cout << "[debug]   probe edge57"
                  << " owner=" << e57.ownerElem
                  << " neigh=" << e57.neighborElem
                  << " bnd=" << (boundary57 ? 1 : 0)
                  << " L=" << e57.edgeLength
                  << " spec=" << e57.spectralRadius
                  << " F3=" << e57.flux[3]
                  << " edgeRE(cell11)=" << re57c11
                  << " edgeRE(cell12)=" << re57c12
                  << "\n";
    } else {
        std::cout << "[debug]   probe edge57 is out of range (edges=" << edges.size() << ").\n";
    }

    std::cout << "[debug] it=" << iteration_
              << " t=" << time_
              << " dt=" << dtUsed
              << " ||R||2=" << normR
              << " rho[min,max]=(" << (*minmaxRho.first)[0] << "," << (*minmaxRho.second)[0] << ")"
              << " p[min,max]=(" << pMin << "@" << pMinCell << "," << pMax << "@" << pMaxCell << ")"
              << " pMin(x,y)=(" << pMinCenter[0] << "," << pMinCenter[1] << ")"
              << " max|flux|=" << fluxMax << "@edge" << fluxMaxEdge
              << " fluxMax(x,y)=(" << fluxMaxCenter[0] << "," << fluxMaxCenter[1] << ")"
              << " max|R|=" << rMax << "@cell" << rMaxCell << ",k" << rMaxComp
              << " max|F_mass,wall|=" << maxWallBoundaryMassFlux()
              << " edges=" << edges.size() << "\n";

    if (dtLocal != nullptr && !dtLocal->empty()) {
        const auto mmDt = std::minmax_element(dtLocal->begin(), dtLocal->end());
        std::cout << "[debug]   dtLocal[min,max]=(" << *mmDt.first << "," << *mmDt.second << ")\n";
    }
}

double FirstorderEuler::maxWallBoundaryMassFlux() const {
    if (boundaryFaces_.empty()) {
        return 0.0;
    }

    const EulerBoundaryConditions bcModel = makeBoundaryConditions(config_);
    const auto flux = makeFlux(config_.fluxScheme);

    double maxAbsMassFlux = 0.0;
    for (const auto& f : boundaryFaces_) {
        const BoundaryKind kind = boundaryKindFromTitle(f.boundaryTitle, config_);
        if (kind != BoundaryKind::WallSlip) {
            continue;
        }

        EulerBoundaryConditions::Context ctx;
        ctx.time = time_;
        ctx.faceCenter = f.center;

        const Conserved UL = U_.at(f.elem);
        const Vec2 nUnit = normalized(f.normal);
        const Eigen::Vector2d n(nUnit[0], nUnit[1]);
        const Conserved UR = bcModel.boundaryState(toBoundaryType(kind), UL, nUnit, ctx);
        const Conserved F = fromEigen((*flux)(toEigen(UL), toEigen(UR), config_.gamma, n));
        maxAbsMassFlux = std::max(maxAbsMassFlux, std::abs(F[0]));
    }
    return maxAbsMassFlux;
}

double FirstorderEuler::l2Norm(const std::vector<Conserved>& values) {
    double acc = 0.0;
    for (const auto& vi : values) {
        for (double x : vi) {
            if (!std::isfinite(x)) {
                return std::numeric_limits<double>::infinity();
            }
            acc += x * x;
        }
    }
    return std::sqrt(acc);
}
