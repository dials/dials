#include <dx2/crystal.h>
#include <gtest/gtest.h>

using Eigen::Matrix3d;
using Eigen::Vector3d;

// A very basic test that we can make a crystal, tests tht gemmi & Eigen is
// installed correctly too.
TEST(ExampleTests, CrystalTest) {
  Vector3d a = {1.0, 0.0, 0.0};
  Vector3d b = {0.0, 1.0, 0.0};
  Vector3d c = {0.0, 0.0, 2.0};
  gemmi::SpaceGroup space_group = *gemmi::find_spacegroup_by_name("P1");
  Crystal crystal(a, b, c, space_group);
  Matrix3d expected_A{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 0.5}};
  EXPECT_EQ(crystal.get_A_matrix(), expected_A);
}
