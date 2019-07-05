#include <scitbx/constants.h>
#include <scitbx/math/r3_rotation.h>
#include <dials/algorithms/spot_prediction/reeke_index_generator.h>
#include <scitbx/array_family/simple_io.h>

using namespace dials::algorithms::reeke_detail;
using namespace dials::algorithms;

using scitbx::mat3;
using scitbx::vec2;
using scitbx::vec3;
using scitbx::math::r3_rotation::axis_and_angle_as_matrix;

template <typename T>
bool almost_equal(T a, T b, T eps) {
  return std::abs(a - b) < eps;
}

template <typename T>
bool almost_equal(vec2<T> a, vec2<T> b, T eps) {
  return std::abs(a[0] - b[0]) < eps && std::abs(a[1] - b[1]);
}

class ReekeInput {
public:
  ReekeInput() {
    dmin = 1.5;

    // Axis and source vectors
    axis = vec3<double>(0, 1, 0);
    source = vec3<double>(0, 0, -1);

    double a = 50;
    ub_beg = mat3<double>(1.0 / a, 0.0, 0.0, 0.0, 1.0 / a, 0.0, 0.0, 0.0, 1.0 / a);

    mat3<double> r_osc = axis_and_angle_as_matrix(axis, 1.0, true);
    ub_end = r_osc * ub_beg;

    margin = 1;
  }

  double dmin;
  vec3<double> axis;
  vec3<double> source;
  mat3<double> ub_beg;
  mat3<double> ub_end;
  std::size_t margin;
};

class TestReekePermutation {
public:
  static void run() {
    ReekeInput input;

    permute_matrix perm(input.ub_beg, input.ub_end, input.axis, input.source);

    // Assert the indices are the right order
    DIALS_ASSERT(perm.index[0] == 2);
    DIALS_ASSERT(perm.index[1] == 0);
    DIALS_ASSERT(perm.index[2] == 1);

    // Assert the permutation is correct
    DIALS_ASSERT(perm.permutation[0] == 0);
    DIALS_ASSERT(perm.permutation[1] == 1);
    DIALS_ASSERT(perm.permutation[2] == 0);
    DIALS_ASSERT(perm.permutation[3] == 0);
    DIALS_ASSERT(perm.permutation[4] == 0);
    DIALS_ASSERT(perm.permutation[5] == 1);
    DIALS_ASSERT(perm.permutation[6] == 1);
    DIALS_ASSERT(perm.permutation[7] == 0);
    DIALS_ASSERT(perm.permutation[8] == 0);

    // Assert the rlv beg and end are ok
    DIALS_ASSERT(almost_equal(perm.rlv_beg[0], 0.0, 1e-7));
    DIALS_ASSERT(almost_equal(perm.rlv_beg[1], 0.0, 1e-7));
    DIALS_ASSERT(almost_equal(perm.rlv_beg[2], 0.02, 1e-7));
    DIALS_ASSERT(almost_equal(perm.rlv_beg[3], 0.02, 1e-7));
    DIALS_ASSERT(almost_equal(perm.rlv_beg[4], 0.0, 1e-7));
    DIALS_ASSERT(almost_equal(perm.rlv_beg[5], 0.0, 1e-7));
    DIALS_ASSERT(almost_equal(perm.rlv_beg[6], 0.0, 1e-7));
    DIALS_ASSERT(almost_equal(perm.rlv_beg[7], 0.02, 1e-7));
    DIALS_ASSERT(almost_equal(perm.rlv_beg[8], 0.0, 1e-7));

    DIALS_ASSERT(almost_equal(perm.rlv_end[0], 0.00034904812874567024, 1e-7));
    DIALS_ASSERT(almost_equal(perm.rlv_end[1], 0.0, 1e-7));
    DIALS_ASSERT(almost_equal(perm.rlv_end[2], 0.019996953903127827, 1e-7));
    DIALS_ASSERT(almost_equal(perm.rlv_end[3], 0.019996953903127827, 1e-7));
    DIALS_ASSERT(almost_equal(perm.rlv_end[4], 0.0, 1e-7));
    DIALS_ASSERT(almost_equal(perm.rlv_end[5], -0.00034904812874567024, 1e-7));
    DIALS_ASSERT(almost_equal(perm.rlv_end[6], 0.0, 1e-7));
    DIALS_ASSERT(almost_equal(perm.rlv_end[7], 0.02, 1e-7));
    DIALS_ASSERT(almost_equal(perm.rlv_end[8], 0.0, 1e-7));

    std::cout << "OK" << std::endl;
  }
};

class TestReekeWithConstantP {
public:
  static void run() {
    ReekeInput input;

    permute_matrix perm(input.ub_beg, input.ub_end, input.axis, input.source);
    compute_constant_with_p cp(
      perm.rlv_beg, perm.rlv_end, input.axis, input.source, input.source);

    DIALS_ASSERT(almost_equal(cp.cp_beg[0], 0.0004, 1e-7));
    DIALS_ASSERT(almost_equal(cp.cp_beg[1], 0.0, 1e-7));
    DIALS_ASSERT(almost_equal(cp.cp_beg[2], -8e-6, 1e-7));
    DIALS_ASSERT(almost_equal(cp.cp_beg[3], -1.6e-7, 1e-7));
    DIALS_ASSERT(almost_equal(cp.cp_beg[4], 0.0, 1e-7));
    DIALS_ASSERT(almost_equal(cp.cp_beg[5], 0.0, 1e-7));
    DIALS_ASSERT(almost_equal(cp.cp_beg[6], -1.6e-7, 1e-7));
    DIALS_ASSERT(almost_equal(cp.cp_beg[7], 0.0, 1e-7));
    DIALS_ASSERT(almost_equal(cp.cp_beg[8], 0.0, 1e-7));
    DIALS_ASSERT(almost_equal(cp.cp_beg[9], 0.0004, 1e-7));
    DIALS_ASSERT(almost_equal(cp.cp_beg[10], 0.0004, 1e-7));
    DIALS_ASSERT(almost_equal(cp.cp_beg[11], 0.0, 1e-7));
    DIALS_ASSERT(almost_equal(cp.cp_beg[12], 0.0, 1e-7));
    DIALS_ASSERT(almost_equal(cp.cp_beg[13], 0.0, 1e-7));
    DIALS_ASSERT(almost_equal(cp.cp_beg[14], 0.04, 1e-7));

    std::cout << "OK" << std::endl;
  }
};

class TestReekeModel {
public:
  static void run() {
    ReekeInput input;

    ReekeModel model(input.ub_beg,
                     input.ub_end,
                     input.axis,
                     input.source,
                     input.source,
                     input.dmin,
                     input.margin);

    DIALS_ASSERT(almost_equal(model.ewald_sphere_p_limits().first[0], -100.0, 1e-7));
    DIALS_ASSERT(almost_equal(model.ewald_sphere_p_limits().first[1], 0.0, 1e-7));
    DIALS_ASSERT(
      almost_equal(model.ewald_sphere_p_limits().second[0], -99.99238475781956, 1e-7));
    DIALS_ASSERT(almost_equal(
      model.ewald_sphere_p_limits().second[1], 0.007615242180436521, 1e-7));
    DIALS_ASSERT(
      almost_equal(model.resolution_p_limits().first[0], -11.11111111111111, 1e-7));
    DIALS_ASSERT(
      almost_equal(model.resolution_p_limits().first[1], -11.11111111111111, 1e-7));
    DIALS_ASSERT(
      almost_equal(model.resolution_p_limits().second[0], -11.657895054618812, 1e-7));
    DIALS_ASSERT(
      almost_equal(model.resolution_p_limits().second[1], -10.5609426155232, 1e-7));

    vec2<int> p_lim = model.p_limits();
    DIALS_ASSERT(p_lim[0] == -12);
    DIALS_ASSERT(p_lim[1] == 2);

    std::cout << "OK" << std::endl;
    DIALS_ASSERT(model.q_limits(-12) == vec2<int>(-32, 33));
    DIALS_ASSERT(model.q_limits(-11) == vec2<int>(-32, 33));
    DIALS_ASSERT(model.q_limits(-10) == vec2<int>(-31, 32));
    DIALS_ASSERT(model.q_limits(-9) == vec2<int>(-29, 31));
    DIALS_ASSERT(model.q_limits(-8) == vec2<int>(-28, 30));
    DIALS_ASSERT(model.q_limits(-7) == vec2<int>(-26, 28));
    DIALS_ASSERT(model.q_limits(-6) == vec2<int>(-24, 26));
    DIALS_ASSERT(model.q_limits(-5) == vec2<int>(-22, 24));
    DIALS_ASSERT(model.q_limits(-4) == vec2<int>(-20, 22));
    DIALS_ASSERT(model.q_limits(-3) == vec2<int>(-18, 19));
    DIALS_ASSERT(model.q_limits(-2) == vec2<int>(-15, 16));
    DIALS_ASSERT(model.q_limits(-1) == vec2<int>(-10, 12));
    DIALS_ASSERT(model.q_limits(0) == vec2<int>(-1, 3));

    std::cout << "OK" << std::endl;

    DIALS_ASSERT(model.r_limits(-8, 19)[0] == vec2<int>(-21, -17));
    DIALS_ASSERT(model.r_limits(-8, 19)[1] == vec2<int>(18, 22));
    DIALS_ASSERT(model.r_limits(-8, 20)[0] == vec2<int>(-20, -16));
    DIALS_ASSERT(model.r_limits(-8, 20)[1] == vec2<int>(17, 21));
    DIALS_ASSERT(model.r_limits(-8, 21)[0] == vec2<int>(-19, -15));
    DIALS_ASSERT(model.r_limits(-8, 21)[1] == vec2<int>(16, 20));
    DIALS_ASSERT(model.r_limits(-8, 22)[0] == vec2<int>(-18, -13));
    DIALS_ASSERT(model.r_limits(-8, 22)[1] == vec2<int>(14, 19));
    DIALS_ASSERT(model.r_limits(-8, 23)[0] == vec2<int>(-16, -12));
    DIALS_ASSERT(model.r_limits(-8, 23)[1] == vec2<int>(13, 17));
    DIALS_ASSERT(model.r_limits(-8, 24)[0] == vec2<int>(-15, -10));
    DIALS_ASSERT(model.r_limits(-8, 24)[1] == vec2<int>(11, 16));
    DIALS_ASSERT(model.r_limits(-8, 25)[0] == vec2<int>(-13, -8));
    DIALS_ASSERT(model.r_limits(-8, 25)[1] == vec2<int>(9, 14));
    DIALS_ASSERT(model.r_limits(-8, 26)[0] == vec2<int>(-11, -5));
    DIALS_ASSERT(model.r_limits(-8, 26)[1] == vec2<int>(6, 12));
    DIALS_ASSERT(model.r_limits(-8, 27)[0] == vec2<int>(-8, 0));
    DIALS_ASSERT(model.r_limits(-8, 27)[1] == vec2<int>(1, 9));
    DIALS_ASSERT(model.r_limits(-8, 28)[0] == vec2<int>(-1, 2));
    DIALS_ASSERT(model.r_limits(-8, 28).size() == 1);
    DIALS_ASSERT(model.r_limits(-7, -25)[0] == vec2<int>(-6, 7));
    DIALS_ASSERT(model.r_limits(-7, -25).size() == 1);
    DIALS_ASSERT(model.r_limits(-7, -24)[0] == vec2<int>(-9, -3));
    DIALS_ASSERT(model.r_limits(-7, -24)[1] == vec2<int>(4, 10));
    DIALS_ASSERT(model.r_limits(-7, -23)[0] == vec2<int>(-12, -7));
    DIALS_ASSERT(model.r_limits(-7, -23)[1] == vec2<int>(8, 13));
    DIALS_ASSERT(model.r_limits(-7, -22)[0] == vec2<int>(-13, -9));
    DIALS_ASSERT(model.r_limits(-7, -22)[1] == vec2<int>(10, 14));
    DIALS_ASSERT(model.r_limits(-7, -21)[0] == vec2<int>(-15, -11));
    DIALS_ASSERT(model.r_limits(-7, -21)[1] == vec2<int>(12, 16));
    DIALS_ASSERT(model.r_limits(-7, -20)[0] == vec2<int>(-16, -12));
    DIALS_ASSERT(model.r_limits(-7, -20)[1] == vec2<int>(13, 17));

    std::cout << "OK" << std::endl;
  }
};

int main() {
  TestReekePermutation::run();
  TestReekeWithConstantP::run();
  TestReekeModel::run();

  return 0;
}
