#include "LinearElement.h"

LinearElement::LinearElement(const Eigen::Vector3i& pointID, const Eigen::Vector3i& faceID, double area, const Eigen::Matrix2d& J):
    Element(pointID, faceID, area),
    _J(J),
    _detJ(J.determinant())
{}