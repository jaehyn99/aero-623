#include "CurvedElement.h"

CurvedElement::CurvedElement(const Eigen::Vector3i& pointID, const Eigen::Vector3i& faceID, double area):
    Element(pointID, faceID, area),
    _internal{0.0, 0.0}
{}