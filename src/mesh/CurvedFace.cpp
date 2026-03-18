#include "CurvedFace.h"

CurvedFace::CurvedFace(const Eigen::Vector2i& pointID, double length, int Q, std::string title):
    Face(pointID, length, title)
{}