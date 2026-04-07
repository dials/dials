// detector_derivatives.h
// Fused C++ kernel for _detector_derivatives (design doc B1).
// Replaces the per-parameter flex arithmetic in
// PredictionParameterisation._detector_derivatives.
//
// The kernel fuses the three-step chain:
//   T1_i = der_i * (-1.0)          // mat3 scalar mul
//   T2_i = D_i   * T1_i            // mat3 * mat3
//   dpv_i = T2_i * pv_i            // mat3 * vec3
// into a single pass over the reflection arrays, eliminating two
// intermediate flex allocations per parameter.
//
// Additionally fuses the subsequent dX/dY computation:
//   dX_i = w_inv_i * (dpv_i.x - dpv_i.z * u_w_inv_i)
//   dY_i = w_inv_i * (dpv_i.y - dpv_i.z * v_w_inv_i)
//
// Operation order mirrors the Python flex left-to-right associativity
// exactly.  See fused_detector_derivatives_design.md for bit-identity
// analysis.

#ifndef DIALS_REFINEMENT_DETECTOR_DERIVATIVES_H
#define DIALS_REFINEMENT_DETECTOR_DERIVATIVES_H

#include <boost/python.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/list.hpp>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/error.h>
#include <scitbx/array_family/ref.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/mat3.h>
#include <scitbx/vec3.h>

namespace dials { namespace refinement {

  using scitbx::mat3;
  using scitbx::vec3;
  namespace af = scitbx::af;

  // Result for one detector parameter: dX/dp and dY/dp over n reflections.
  // Empty dX/dY (size 0) indicates a null (None) parameter.
  struct dXY_pair {
    af::shared<double> dX;
    af::shared<double> dY;
  };

  // Fused detector-derivative kernel.
  //
  // For each non-null parameter p and each reflection i:
  //   neg_der = -der_p_i                          // mat3 sign flip
  //   M       = D_i * neg_der                     // mat3 * mat3
  //   dpv     = M * pv_i                          // mat3 * vec3
  //   dX_i    = w_inv_i * (dpv.x - dpv.z * u_w_inv_i)
  //   dY_i    = w_inv_i * (dpv.y - dpv.z * v_w_inv_i)
  //
  // dd_ddet_p is a Python list of length nparam.  Each entry is either:
  //   - None      => null parameter; result entry gets empty dX/dY arrays.
  //   - flex.mat3_double of length n (dense, broadcast) => active parameter.
  //
  // D, pv, w_inv, u_w_inv, v_w_inv are all length n (sub-panel selection).
  //
  // Returns af::shared<dXY_pair> of length nparam.
  //
  // Bit-identity: arithmetic order exactly mirrors the Python chain
  //   (D * (der * -1.0)) * pv => (D * neg_der) * pv.
  // No reassociation is performed.
  inline af::shared<dXY_pair> compute_detector_derivatives(
    af::const_ref<mat3<double> > const& D,
    af::const_ref<vec3<double> > const& pv,
    af::const_ref<double> const& w_inv,
    af::const_ref<double> const& u_w_inv,
    af::const_ref<double> const& v_w_inv,
    boost::python::object const& dd_ddet_p)  // list of flex.mat3_double or None
  {
    namespace bp = boost::python;

    std::size_t n = D.size();
    DIALS_ASSERT(pv.size() == n);
    DIALS_ASSERT(w_inv.size() == n);
    DIALS_ASSERT(u_w_inv.size() == n);
    DIALS_ASSERT(v_w_inv.size() == n);

    std::size_t nparam = bp::len(dd_ddet_p);
    af::shared<dXY_pair> result(nparam);

    for (std::size_t p = 0; p < nparam; ++p) {
      bp::object entry = dd_ddet_p[p];

      if (entry.ptr() == Py_None) {
        // Null parameter: leave dX/dY empty (size 0).
        continue;
      }

      // Extract dense flex.mat3_double for this parameter.
      af::const_ref<mat3<double> > der =
        bp::extract<af::const_ref<mat3<double> > >(entry);
      DIALS_ASSERT(der.size() == n);

      af::shared<double> dX(n);
      af::shared<double> dY(n);

      for (std::size_t i = 0; i < n; ++i) {
        // Replicate Python:  (D * (der * -1.0)) * pv
        // Step 1: neg_der = der_i * (-1.0)  [mat3 scalar mul]
        // Step 2: M       = D_i   * neg_der [mat3 * mat3]
        // Step 3: dpv     = M     * pv_i    [mat3 * vec3]
        // Preserved left-to-right: (D * (der * -1.0)) * pv
        // which equals D*(neg_der)*pv, not -(D*der*pv).  Same bits.
        mat3<double> neg_der = der[i] * (-1.0);
        mat3<double> M = D[i] * neg_der;
        vec3<double> dpv = M * pv[i];

        double du = dpv[0];
        double dv = dpv[1];
        double dw = dpv[2];

        // Quotient-rule expansion — preserve operand order exactly.
        dX[i] = w_inv[i] * (du - dw * u_w_inv[i]);
        dY[i] = w_inv[i] * (dv - dw * v_w_inv[i]);
      }

      result[p].dX = dX;
      result[p].dY = dY;
    }

    return result;
  }

}}  // namespace dials::refinement
#endif  // DIALS_REFINEMENT_DETECTOR_DERIVATIVES_H
