// xl_derivatives.h
// Fused C++ kernel for _xl_derivatives (design doc opp 2.3).
// Replaces the per-parameter inner loop in
// XYPhiPredictionParameterisation._xl_derivatives.
//
// Operation order mirrors the Python flex left-to-right associativity exactly.
// See fused_xl_derivatives_design.md for the full bit-identity analysis.

#ifndef DIALS_REFINEMENT_XL_DERIVATIVES_H
#define DIALS_REFINEMENT_XL_DERIVATIVES_H

#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/ref.h>
#include <scitbx/vec3.h>
#include <scitbx/mat3.h>
#include <dials/error.h>

namespace dials { namespace refinement {

  using scitbx::mat3;
  using scitbx::vec3;
  namespace af = scitbx::af;

  // Output struct: one parameter's contribution over N reflections.
  struct XlDerivResult {
    af::shared<vec3<double> > dpv;  // size N
    af::shared<double> dphi;        // size N
  };

  // Fused computation for a SINGLE parameter's derivative matrix `der`.
  //
  // When b_matrix == true:
  //   tmp_i = fixed_rotation_i * ((der * B_or_U_i) * h_i)   [left-to-right]
  // When b_matrix == false:
  //   tmp_i = fixed_rotation_i * ((B_or_U_i * der) * h_i)   [left-to-right]
  //
  // Then for each reflection i:
  //   rotated_i = tmp_i.unit_rotate_around_origin(axis_i.normalize(), phi_calc_i)
  //   dr_i      = setting_rotation_i * rotated_i
  //   dphi_i    = (-1.0 * (dr_i * s1_i)) / e_r_s0_i
  //   dpv_i     = D_i * (dr_i + e_X_r_i * dphi_i)
  //
  // All per-reflection arrays must have the same length N.
  // Matrix products are performed strictly left-to-right to preserve
  // bit-identity with the Python/flex path.
  inline XlDerivResult compute_xl_derivative_one(
    mat3<double> const& der,
    af::const_ref<mat3<double> > const& fixed_rotation,
    af::const_ref<mat3<double> > const& setting_rotation,
    af::const_ref<mat3<double> > const& B_or_U,  // B if b_matrix, else U
    af::const_ref<vec3<double> > const& h,
    af::const_ref<vec3<double> > const& axis,
    af::const_ref<double> const& phi_calc,
    af::const_ref<vec3<double> > const& s1,
    af::const_ref<vec3<double> > const& e_X_r,
    af::const_ref<double> const& e_r_s0,
    af::const_ref<mat3<double> > const& D,
    bool b_matrix) {
    std::size_t N = h.size();
    DIALS_ASSERT(fixed_rotation.size() == N);
    DIALS_ASSERT(setting_rotation.size() == N);
    DIALS_ASSERT(B_or_U.size() == N);
    DIALS_ASSERT(axis.size() == N);
    DIALS_ASSERT(phi_calc.size() == N);
    DIALS_ASSERT(s1.size() == N);
    DIALS_ASSERT(e_X_r.size() == N);
    DIALS_ASSERT(e_r_s0.size() == N);
    DIALS_ASSERT(D.size() == N);

    XlDerivResult out;
    out.dpv.resize(N);
    out.dphi.resize(N);

    for (std::size_t i = 0; i < N; ++i) {
      // Step 1: tmp = fixed_rotation_i * (der op B_or_U_i) * h_i [left-to-right]
      vec3<double> tmp;
      if (b_matrix) {
        // Python: der * B * h  =>  (der * B_i) * h_i
        mat3<double> dB = der * B_or_U[i];  // der * B_i
        vec3<double> v = dB * h[i];         // (der * B_i) * h_i
        tmp = fixed_rotation[i] * v;
      } else {
        // Python: U * der * h  =>  (U_i * der) * h_i
        mat3<double> Ud = B_or_U[i] * der;  // U_i * der
        vec3<double> v = Ud * h[i];         // (U_i * der) * h_i
        tmp = fixed_rotation[i] * v;
      }

      // Step 2: rotate. Mirror flex_vec3_double.cpp:130-145 exactly:
      //   normalize axis per-element, then unit_rotate_around_origin.
      // Do NOT hoist normalize() out of the loop — doing so changes rounding.
      vec3<double> unit = axis[i].normalize();
      vec3<double> rotated = tmp.unit_rotate_around_origin(unit, phi_calc[i]);

      // Step 3: dr = setting_rotation_i * rotated
      vec3<double> dr = setting_rotation[i] * rotated;

      // Step 4: dphi = (-1.0 * (dr * s1_i)) / e_r_s0_i
      // Python: -1.0 * dr.dot(s1) / e_r_s0  =>  ((-1.0 * dot) / e_r_s0)
      // scitbx vec3: operator*(vec3, vec3) is the dot product (vec3.h:413).
      // Left-to-right: multiply by -1.0 first, THEN divide. Do not negate after
      // dividing (would change associativity and break bit-identity).
      double dot = dr * s1[i];
      double dphi_i = (-1.0 * dot) / e_r_s0[i];

      // Step 5: dpv = D_i * (dr + e_X_r_i * dphi_i)
      vec3<double> e_X_r_dphi = e_X_r[i] * dphi_i;
      vec3<double> sum = dr + e_X_r_dphi;
      vec3<double> dpv_i = D[i] * sum;

      out.dphi[i] = dphi_i;
      out.dpv[i] = dpv_i;
    }
    return out;
  }

}}  // namespace dials::refinement
#endif  // DIALS_REFINEMENT_XL_DERIVATIVES_H
