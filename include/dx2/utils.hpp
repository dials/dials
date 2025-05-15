#ifndef DX2_UTILS_H
#define DX2_UTILS_H
#include <Eigen/Dense>
#include <math.h>

using Eigen::Vector3d;

double angle_between_vectors_degrees(Vector3d v1, Vector3d v2) {
  double l1 = v1.norm();
  double l2 = v2.norm();
  double dot = v1.dot(v2);
  double normdot = dot / (l1 * l2);
  if (std::abs(normdot - 1.0) < 1E-6) {
    return 0.0;
  }
  if (std::abs(normdot + 1.0) < 1E-6) {
    return 180.0;
  }
  double angle = std::acos(normdot) * 180.0 / M_PI;
  return angle;
}

#endif // DX2_UTILS_H
