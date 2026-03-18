#include "Element.h"

Element::Element(const Eigen::Vector3i& pointID, const Eigen::Vector3i& faceID, double area)://, const Eigen::Vector2d& centroid):
    _pointID(pointID),
    _faceID(faceID),
    _area(area)
    //_centroid(centroid)
{}