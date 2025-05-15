#include <dx2/crystal.hpp>
#include <gtest/gtest.h>

using Eigen::Matrix3d;
using Eigen::Vector3d;

// Test that we can correctly construct a crystal
// A useful test that gemmi & Eigen is installed correctly too.
TEST(ModelTests, CrystalTest) {
  // Initialise a crystal
  Vector3d a = {10.0, 1.0, 2.0};
  Vector3d b = {1.0, 20.0, 4.0};
  Vector3d c = {2.0, 2.0, 5.0};
  gemmi::SpaceGroup space_group = *gemmi::find_spacegroup_by_name("P1");
  Crystal crystal(a, b, c, space_group);

  // Expected values based on equivalent dials/dxtbx code.
  std::array<double, 6> expected_cell = {10.246950765959598, 20.42057785666214,
                                         5.744562646538029,  58.09405586165891,
                                         57.06933979197713,  79.53690718052404};
  Matrix3d expected_A{
      {0.10861865407319952, -0.001180637544273908, -0.042502951593860694},
      {0.0035419126328217194, 0.054309327036599776, -0.044864226682408505},
      {-0.0448642266824085, -0.021251475796930347, 0.2349468713105077}};
  Matrix3d expected_B{
      {0.09759000729485333, 0.0, 0.0},
      {-0.018022224347326803, 0.04979825148603455, 0.0},
      {-0.06304558588313347, -0.030375257705251953, 0.24293894720401843}};
  Matrix3d expected_U{
      {0.9759000729485331, -0.13042399198723334, -0.1749532221285499},
      {0.09759000729485345, 0.9779428053733645, -0.184672845580136},
      {0.19518001458970666, 0.16314855724948435, 0.9671025334328175}};

  // Now get the matrices and test against expected values
  Matrix3d Amatrix = crystal.get_A_matrix();
  Matrix3d Bmatrix = crystal.get_B_matrix();
  Matrix3d Umatrix = crystal.get_U_matrix();
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_NEAR(Amatrix(i, j), expected_A(i, j), 1E-12);
    }
  }
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_NEAR(Bmatrix(i, j), expected_B(i, j), 1E-12);
    }
  }
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_NEAR(Umatrix(i, j), expected_U(i, j), 1E-12);
    }
  }
  // Check the cell against expected values
  gemmi::UnitCell cell = crystal.get_unit_cell();
  std::array<double, 6> params = {cell.a,     cell.b,    cell.c,
                                  cell.alpha, cell.beta, cell.gamma};
  for (int i = 0; i < 6; i++) {
    EXPECT_NEAR(params[i], expected_cell[i], 1E-12);
  }

  // Now do niggli reduction, and check the quantities against dials/dxtbx
  // equivalent
  crystal.niggli_reduce();

  gemmi::UnitCell reduced_cell = crystal.get_unit_cell();
  std::array<double, 6> expected_reduced_cell = {
      5.744562646538028, 8.602325267042627, 17.34935157289747,
      98.47679564504885, 92.30016540862056, 91.15952321378222};
  std::array<double, 6> reduced_params = {
      reduced_cell.a,     reduced_cell.b,    reduced_cell.c,
      reduced_cell.alpha, reduced_cell.beta, reduced_cell.gamma};
  for (int i = 0; i < 6; i++) {
    EXPECT_NEAR(reduced_params[i], expected_reduced_cell[i], 1E-12);
  }

  Matrix3d expected_A_reduced{
      {0.06375442739079101, 0.10861865407319951, -0.001180637544273908},
      {0.06729634002361276, 0.003541912632821717, 0.054309327036599776},
      {0.1475796930342385, -0.0448642266824085, -0.021251475796930347}};
  Matrix3d expected_B_reduced{
      {0.17407765595569788, 0.0, 0.0},
      {0.0035233772055836498, 0.11627144778426224, 0.0},
      {0.00759905906137586, 0.01744438086416477, 0.05833114204030993}};
  Matrix3d expected_U_reduced{
      {0.3481553119113958, 0.9372183366852653, -0.020240261084859688},
      {0.3481553119113958, -0.10922469337309479, 0.9310520099035458},
      {0.870388279778489, -0.33119745732486827, -0.3643246995274744}};
  Matrix3d Amatrix_reduced = crystal.get_A_matrix();
  Matrix3d Bmatrix_reduced = crystal.get_B_matrix();
  Matrix3d Umatrix_reduced = crystal.get_U_matrix();
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_NEAR(Amatrix_reduced(i, j), expected_A_reduced(i, j), 1E-12);
    }
  }
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_NEAR(Bmatrix_reduced(i, j), expected_B_reduced(i, j), 1E-12);
    }
  }
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_NEAR(Umatrix_reduced(i, j), expected_U_reduced(i, j), 1E-12);
    }
  }
}
