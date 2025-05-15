#ifndef DX2_MODEL_GONIO_H
#define DX2_MODEL_GONIO_H
#include <Eigen/Dense>
#include <math.h>
#include <nlohmann/json.hpp>

using json = nlohmann::json;
using Eigen::Matrix3d;
using Eigen::Vector3d;

Matrix3d axis_and_angle_as_matrix(Vector3d axis, double angle,
                                  bool deg = false) {
  double q0 = 0.0;
  double q1 = 0.0;
  double q2 = 0.0;
  double q3 = 0.0;
  if (deg) {
    angle *= M_PI / 180.0;
  }
  if (!(std::fmod(angle, 2.0 * M_PI))) {
    q0 = 1.0;
  } else {
    double h = 0.5 * angle;
    q0 = std::cos(h);
    double s = std::sin(h);
    Vector3d n = axis / axis.norm();
    q1 = n[0] * s;
    q2 = n[1] * s;
    q3 = n[2] * s;
  }
  Matrix3d m{{2 * (q0 * q0 + q1 * q1) - 1, 2 * (q1 * q2 - q0 * q3),
              2 * (q1 * q3 + q0 * q2)},
             {2 * (q1 * q2 + q0 * q3), 2 * (q0 * q0 + q2 * q2) - 1,
              2 * (q2 * q3 - q0 * q1)},
             {2 * (q1 * q3 - q0 * q2), 2 * (q2 * q3 + q0 * q1),
              2 * (q0 * q0 + q3 * q3) - 1}};
  return m;
}

class Goniometer {
  // A class to represent a multi-axis goniometer.
public:
  Goniometer() = default;
  Goniometer(std::vector<Vector3d> axes, std::vector<double> angles,
             std::vector<std::string> names, std::size_t scan_axis);
  Goniometer(json goniometer_data);
  Matrix3d get_setting_rotation() const;
  Matrix3d get_sample_rotation() const;
  Vector3d get_rotation_axis() const;
  json to_json() const;

protected:
  void init(); // Sets the matrices from the axes and angles
  Matrix3d calculate_setting_rotation();
  Matrix3d calculate_sample_rotation();
  Matrix3d sample_rotation_{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};  // F
  Vector3d rotation_axis_{{1.0, 0.0, 0.0}};                    // R'
  Matrix3d setting_rotation_{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}; // S
  std::vector<Vector3d> axes_{{1.0, 0.0, 0.0}};
  std::vector<double> angles_{{0.0}};
  std::vector<std::string> names_{"omega"};
  std::size_t scan_axis_{0};
};

void Goniometer::init() {
  // Sets the matrices from the axes and angles
  setting_rotation_ = calculate_setting_rotation();
  sample_rotation_ = calculate_sample_rotation();
  rotation_axis_ = axes_[scan_axis_];
}

Matrix3d Goniometer::calculate_setting_rotation() {
  Matrix3d setting_rotation{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  for (std::size_t i = scan_axis_ + 1; i < axes_.size(); i++) {
    Matrix3d R = axis_and_angle_as_matrix(axes_[i], angles_[i], true);
    setting_rotation = R * setting_rotation;
  }
  return setting_rotation;
}

Matrix3d Goniometer::calculate_sample_rotation() {
  Matrix3d sample_rotation{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  for (std::size_t i = 0; i < scan_axis_; i++) {
    Matrix3d R = axis_and_angle_as_matrix(axes_[i], angles_[i], true);
    sample_rotation = R * sample_rotation;
  }
  return sample_rotation;
}

Goniometer::Goniometer(std::vector<Vector3d> axes, std::vector<double> angles,
                       std::vector<std::string> names, std::size_t scan_axis)
    : axes_{axes.begin(), axes.end()}, angles_{angles.begin(), angles.end()},
      names_{names.begin(), names.end()}, scan_axis_{scan_axis} {
  if (scan_axis_ >= axes_.size()) {
    throw std::invalid_argument(
        "Goniometer scan axis number is out of range of axis length");
  }
  init();
}

Matrix3d Goniometer::get_setting_rotation() const { return setting_rotation_; }

Matrix3d Goniometer::get_sample_rotation() const { return sample_rotation_; }

Vector3d Goniometer::get_rotation_axis() const { return rotation_axis_; }

Goniometer::Goniometer(json goniometer_data) {
  std::vector<std::string> required_keys = {"axes", "angles", "names",
                                            "scan_axis"};
  for (const auto &key : required_keys) {
    if (goniometer_data.find(key) == goniometer_data.end()) {
      throw std::invalid_argument("Key " + key +
                                  " is missing from the input goniometer JSON");
    }
  }
  std::vector<Vector3d> axes;
  for (json::iterator it = goniometer_data["axes"].begin();
       it != goniometer_data["axes"].end(); it++) {
    Vector3d axis;
    axis[0] = (*it)[0];
    axis[1] = (*it)[1];
    axis[2] = (*it)[2];
    axes.push_back(axis);
  }
  axes_ = axes;
  std::vector<double> angles;
  for (json::iterator it = goniometer_data["angles"].begin();
       it != goniometer_data["angles"].end(); it++) {
    angles.push_back(*it);
  }
  angles_ = angles;
  std::vector<std::string> names;
  for (json::iterator it = goniometer_data["names"].begin();
       it != goniometer_data["names"].end(); it++) {
    names.push_back(*it);
  }
  names_ = names;
  scan_axis_ = goniometer_data["scan_axis"];
  init();
}

json Goniometer::to_json() const {
  json goniometer_data;
  goniometer_data["axes"] = axes_;
  goniometer_data["angles"] = angles_;
  goniometer_data["names"] = names_;
  goniometer_data["scan_axis"] = scan_axis_;
  return goniometer_data;
}

#endif // DX2_MODEL_GONIO_H
