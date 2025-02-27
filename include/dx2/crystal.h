#ifndef DX2_MODEL_CRYSTAL_H
#define DX2_MODEL_CRYSTAL_H
#include "utils.h"
#include <Eigen/Dense>
#include <gemmi/cellred.hpp>
#include <gemmi/math.hpp>
#include <gemmi/symmetry.hpp>
#include <gemmi/unitcell.hpp>
#include <iostream>
#include <math.h>
#include <nlohmann/json.hpp>

using json = nlohmann::json;
using Eigen::Matrix3d;
using Eigen::Vector3d;

Matrix3d Matrix3d_from_gemmi_cb(gemmi::Op cb) {
  std::array<std::array<int, 3>, 3> rot = cb.rot;
  return Matrix3d{{(double)(rot[0][0] / cb.DEN), (double)(rot[1][0] / cb.DEN),
                   (double)(rot[2][0] / cb.DEN)},
                  {(double)(rot[0][1] / cb.DEN), (double)(rot[1][1] / cb.DEN),
                   (double)(rot[2][1] / cb.DEN)},
                  {(double)(rot[0][2] / cb.DEN), (double)(rot[1][2] / cb.DEN),
                   (double)(rot[2][2] / cb.DEN)}};
}

class Crystal {
public:
  Crystal() = default;
  Crystal(Vector3d a, Vector3d b, Vector3d c, gemmi::SpaceGroup space_group);
  Crystal(json crystal_data);
  gemmi::UnitCell get_unit_cell() const;
  gemmi::SpaceGroup get_space_group() const;
  Matrix3d get_A_matrix() const;
  Matrix3d get_U_matrix() const;
  Matrix3d get_B_matrix() const;
  void niggli_reduce();
  json to_json() const;

protected:
  void init_from_abc(Vector3d a, Vector3d b, Vector3d c);
  gemmi::SpaceGroup space_group_;
  gemmi::UnitCell unit_cell_;
  Matrix3d B_;
  Matrix3d A_;
  Matrix3d U_;
};

void Crystal::init_from_abc(Vector3d a, Vector3d b, Vector3d c) {
  // calculate B matrix, A matrix, set the input cell values
  Matrix3d A{{a[0], a[1], a[2]}, {b[0], b[1], b[2]}, {c[0], c[1], c[2]}};
  A_ = A.inverse();
  double real_space_a = a.norm();
  double real_space_b = b.norm();
  double real_space_c = c.norm();
  double alpha = angle_between_vectors_degrees(b, c);
  double beta = angle_between_vectors_degrees(c, a);
  double gamma = angle_between_vectors_degrees(a, b);
  unit_cell_ = {real_space_a, real_space_b, real_space_c, alpha, beta, gamma};
  gemmi::Mat33 B = unit_cell_.frac.mat;
  B_ << B.a[0][0], B.a[1][0], B.a[2][0], B.a[0][1], B.a[1][1], B.a[2][1],
      B.a[0][2], B.a[1][2],
      B.a[2][2]; // Transpose due to different definition of B form
  U_ = A_ * B_.inverse();
}

Crystal::Crystal(Vector3d a, Vector3d b, Vector3d c,
                 gemmi::SpaceGroup space_group)
    : space_group_(space_group) {
  init_from_abc(a, b, c);
}

Crystal::Crystal(json crystal_data) {
  std::vector<std::string> required_keys = {"real_space_a", "real_space_b",
                                            "real_space_c",
                                            "space_group_hall_symbol"};
  for (const auto &key : required_keys) {
    if (crystal_data.find(key) == crystal_data.end()) {
      throw std::invalid_argument("Key " + key +
                                  " is missing from the input crystal JSON");
    }
  }
  Vector3d rsa{{crystal_data["real_space_a"][0],
                crystal_data["real_space_a"][1],
                crystal_data["real_space_a"][2]}};
  Vector3d rsb{{crystal_data["real_space_b"][0],
                crystal_data["real_space_b"][1],
                crystal_data["real_space_b"][2]}};
  Vector3d rsc{{crystal_data["real_space_c"][0],
                crystal_data["real_space_c"][1],
                crystal_data["real_space_c"][2]}};
  space_group_ =
      *gemmi::find_spacegroup_by_name(crystal_data["space_group_hall_symbol"]);
  init_from_abc(rsa, rsb, rsc);
}

void Crystal::niggli_reduce() {
  char centering{'P'};
  gemmi::GruberVector gv(unit_cell_, centering, true);
  gv.niggli_reduce();
  gemmi::UnitCell niggli_cell = gv.get_cell();
  std::cout << "Input cell:" << std::endl;
  std::cout << unit_cell_.a << " " << unit_cell_.b << " " << unit_cell_.c << " "
            << unit_cell_.alpha << " " << unit_cell_.beta << " "
            << unit_cell_.gamma << std::endl;
  std::cout << "Reduced cell:" << std::endl;
  std::cout << niggli_cell.a << " " << niggli_cell.b << " " << niggli_cell.c
            << " " << niggli_cell.alpha << " " << niggli_cell.beta << " "
            << niggli_cell.gamma << std::endl;
  unit_cell_ = niggli_cell;
  gemmi::Op cb = *gv.change_of_basis;

  Matrix3d cb_op = Matrix3d_from_gemmi_cb(cb);
  A_ = A_ * cb_op.inverse();
  gemmi::Mat33 B = unit_cell_.frac.mat;
  B_ << B.a[0][0], B.a[1][0], B.a[2][0], B.a[0][1], B.a[1][1], B.a[2][1],
      B.a[0][2], B.a[1][2],
      B.a[2][2]; // Transpose due to different definition of B form
  U_ = A_ * B_.inverse();
}

gemmi::UnitCell Crystal::get_unit_cell() const { return unit_cell_; }

gemmi::SpaceGroup Crystal::get_space_group() const { return space_group_; }

Matrix3d Crystal::get_A_matrix() const { return A_; }

Matrix3d Crystal::get_B_matrix() const { return B_; }

Matrix3d Crystal::get_U_matrix() const { return U_; }

json Crystal::to_json() const {
  json crystal_data;
  crystal_data["__id__"] = "crystal";
  Matrix3d A_inv = A_.inverse();
  Vector3d rsa{{A_inv(0, 0), A_inv(0, 1), A_inv(0, 2)}};
  Vector3d rsb{{A_inv(1, 0), A_inv(1, 1), A_inv(1, 2)}};
  Vector3d rsc{{A_inv(2, 0), A_inv(2, 1), A_inv(2, 2)}};
  crystal_data["real_space_a"] = rsa;
  crystal_data["real_space_b"] = rsb;
  crystal_data["real_space_c"] = rsc;
  crystal_data["space_group_hall_symbol"] = space_group_.hall;
  return crystal_data;
}

#endif // DX2_MODEL_CRYSTAL_H
