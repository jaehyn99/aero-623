#include "Face.h"

Face::Face(const Eigen::Vector2i& pointID, double length, std::string title):
    _pointID(pointID),
    _length(length),
    _title(title)
{}

bool Face::operator==(const Face& other) const noexcept{
    return (_pointID[0] == other._pointID[0] && _pointID[1] == other._pointID[1]) || (_pointID[0] == other._pointID[1] && _pointID[1] == other._pointID[0]);
}