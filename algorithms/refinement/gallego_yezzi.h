
#ifndef DIALS_REFINEMENT_GALLEGO_YEZZI_H
#define DIALS_REFINEMENT_GALLEGO_YEZZI_H

#include <cmath>
#include <scitbx/vec3.h>
#include <scitbx/mat3.h>
#include <scitbx/math/r3_rotation.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/error.h>

namespace dials { namespace refinement {

  using scitbx::mat3;
  using scitbx::vec3;
  using scitbx::math::r3_rotation::axis_and_angle_as_matrix;

  mat3<double> skew_symm(vec3<double> v) {
    mat3<double> L1(0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 1.0, 0.0);
    mat3<double> L2(0.0, 0.0, 1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0);
    mat3<double> L3(0.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0);

    return v[0] * L1 + v[1] * L2 + v[2] * L3;
  }

  af::shared<mat3<double> > dRq_de(const af::const_ref<double> &theta,
                                   const af::const_ref<vec3<double> > &e1,
                                   const af::const_ref<vec3<double> > &q) {
    // Calculate the derivative of a rotated vector with respect to the axis
    // of rotation. The result is a 3*3 matrix. This function implements the
    // method of Gallego & Yezzi (equn 8 in http://arxiv.org/pdf/1312.0788.pdf)

    DIALS_ASSERT(theta.size() == e1.size());
    DIALS_ASSERT(theta.size() == q.size());

    af::shared<mat3<double> > result(theta.size(),
                                     af::init_functor_null<mat3<double> >());

    // I(3)
    mat3<double> I3(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);

    // null matrix
    mat3<double> null_mat(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

    for (std::size_t i = 0; i < result.size(); i++) {
      // for angle near zero immediately return null mat
      if (fabs(theta[i]) < 1.e-20) {
        result[i] = null_mat;
        continue;
      }

      // ensure the axis is unit
      vec3<double> e1_u = e1[i].normalize();

      // rotation matrix R
      mat3<double> R = axis_and_angle_as_matrix(e1_u, theta[i]);

      // rotation vector v
      vec3<double> v = theta[i] * e1_u;

      // skew q
      mat3<double> q_x = skew_symm(q[i]);

      // skew v
      mat3<double> v_x = skew_symm(v);

      // outer product, v * v^T
      mat3<double> vvt(v[0] * v[0],
                       v[0] * v[1],
                       v[0] * v[2],
                       v[1] * v[0],
                       v[1] * v[1],
                       v[1] * v[2],
                       v[2] * v[0],
                       v[2] * v[1],
                       v[2] * v[2]);

      // R^T
      mat3<double> Rt = R.transpose();

      // do calculation and put this element in the result
      result[i] = (-1.0 / theta[i]) * R * q_x * (vvt + (Rt - I3) * v_x);
    }

    return result;
  }

}}  // namespace dials::refinement

#endif  // DIALS_REFINEMENT_GALLEGO_YEZZI_H
