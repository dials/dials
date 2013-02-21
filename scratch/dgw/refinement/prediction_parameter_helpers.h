
#ifndef DIALS_REFINEMENT_PREDICTION_PARAMETER_HELPERS_H
#define DIALS_REFINEMENT_PREDICTION_PARAMETER_HELPERS_H

#include <scitbx/vec3.h>
#include <scitbx/mat3.h>
#include <scitbx/array_family/flex_types.h>
#include <cctbx/miller.h>

namespace dials { namespace refinement {

  using scitbx::vec3;
  using scitbx::mat3;
  using scitbx::af::flex_double;

  typedef cctbx::miller::index <> miller_index;
  typedef scitbx::af::flex <vec3 <double> >::type flex_vec3_double;
  typedef scitbx::af::flex <mat3 <double> >::type flex_mat3_double;
  typedef scitbx::af::flex <miller_index >::type flex_miller_index;

  /**
   * Calculate derivatives of pv wrt a single parameter of the FIRST detector
   * parameterisation only.
   */
  vec3 <double> detector_pv_derivative(
      mat3 <double> D, mat3 <double> dd_ddet_p, vec3 <double> pv) {
    return -D * dd_ddet_p * pv;
  }

  /**
   * Calculate derivatives of pv wrt each parameter of the FIRST detector
   * parameterisation only.
   */
  flex_vec3_double detector_pv_derivative(
      mat3 <double> D, const flex_mat3_double &dd_ddet_p, vec3 <double> pv) {
    flex_vec3_double dpv_ddet_p(dd_ddet_p.size());
    for (std::size_t i = 0; i < dd_ddet_p.size(); ++i) {
      dpv_ddet_p[i] = detector_pv_derivative(D, dd_ddet_p[i], pv);
    }
    return dpv_ddet_p;
  }

  /**
   * Calc derivatives of phi wrt each parameter of each source
   * parameterisation that is present.
   */
  double source_phi_derivative(
      vec3 <double> r, vec3 <double> ds0_dsrc_p, double e_r_s0) {
    return -(r * ds0_dsrc_p) / e_r_s0;
  }

  /**
   * Calc derivatives of phi wrt each parameter of each source
   * parameterisation that is present.
   */
  flex_double source_phi_derivative(
      vec3 <double> r, const flex_vec3_double &ds0_dsrc_p, double e_r_s0) {
    flex_double dphi_dsrc_p(ds0_dsrc_p.size());
    for (std::size_t i = 0; i < ds0_dsrc_p.size(); ++i) {
      dphi_dsrc_p[i] = source_phi_derivative(r, ds0_dsrc_p[i], e_r_s0);
    }
    return dphi_dsrc_p;
  }

  /**
   * Calc derivatives of pv wrt each parameter of each source
   * parameterisation that is present.
   */
  vec3 <double> source_pv_derivative(
      mat3 <double> D, vec3 <double> e_X_r, double dphi_dsrc_p, vec3 <double> ds0_dsrc_p) {
    return D * (e_X_r * dphi_dsrc_p + ds0_dsrc_p);
  }

  /**
   * Calc derivatives of pv wrt each parameter of each source
   * parameterisation that is present.
   */
  flex_vec3_double source_pv_derivative(
      mat3 <double> D, vec3 <double> e_X_r, const flex_double &dphi_dsrc_p,
      const flex_vec3_double &ds0_dsrc_p) {
    SCITBX_ASSERT(ds0_dsrc_p.size() == dphi_dsrc_p.size());
    flex_vec3_double dpv_dsrc_p(ds0_dsrc_p.size());
    for (std::size_t i = 0; i < ds0_dsrc_p.size(); ++i) {
      dpv_dsrc_p[i] = source_pv_derivative(
          D, e_X_r, dphi_dsrc_p[i], ds0_dsrc_p[i]);
    }
    return dpv_dsrc_p;
  }

  /**
   * Calc derivatives of r wrt each parameter of each crystal
   * orientation parameterisation that is present.
   */
   vec3 <double> crystal_orientation_r_derivative(
      mat3 <double> R, mat3 <double> dU_dxlo_p, mat3 <double> B,
      miller_index h) {
    return R * dU_dxlo_p * B * h;
  }

  /**
   * Calc derivatives of r wrt each parameter of each crystal
   * orientation parameterisation that is present.
   */
  flex_vec3_double crystal_orientation_r_derivative(
      mat3 <double> R, const flex_mat3_double &dU_dxlo_p, mat3 <double> B,
      miller_index h) {
    flex_vec3_double dr_dxlo_p(dU_dxlo_p.size());
    for (std::size_t i = 0; i < dU_dxlo_p.size(); ++i) {
      dr_dxlo_p[i] = crystal_orientation_r_derivative(
          R, dU_dxlo_p[i], B, h);
    }
    return dr_dxlo_p;
  }

  /**
   * Calc derivatives of phi wrt each parameter of each crystal
   * orientation parameterisation that is present.
   */
  double crystal_orientation_phi_derivative(
      vec3 <double> der, vec3 <double> s, double e_r_s0) {
    return -(der * s) / e_r_s0;
  }

  /**
   * Calc derivatives of phi wrt each parameter of each crystal
   * orientation parameterisation that is present.
   */
  flex_double crystal_orientation_phi_derivative(
      const flex_vec3_double &dr_dxlo_p, vec3 <double> s, double e_r_s0) {
    flex_double dphi_dxlo_p(dr_dxlo_p.size());
    for (std::size_t i = 0; i < dr_dxlo_p.size(); ++i) {
      dphi_dxlo_p[i] = crystal_orientation_phi_derivative(
          dr_dxlo_p[i], s, e_r_s0);
    }
    return dphi_dxlo_p;
  }

  /**
   * Calc derivatives of pv wrt each parameter of each crystal
   * orientation parameterisation that is present.
   */
  vec3 <double> crystal_orientation_pv_derivative(
      mat3 <double> D, vec3 <double> dr_dxlo_p, vec3 <double> e_X_r,
      double dphi_dxlo_p) {
    return D * (dr_dxlo_p + (e_X_r * dphi_dxlo_p));
  }

  /**
   * Calc derivatives of pv wrt each parameter of each crystal
   * orientation parameterisation that is present.
   */
  flex_vec3_double crystal_orientation_pv_derivative(
      mat3 <double> D, const flex_vec3_double &dr_dxlo_p, vec3 <double> e_X_r,
      const flex_double &dphi_dxlo_p) {
    SCITBX_ASSERT(dr_dxlo_p.size() == dphi_dxlo_p.size());
    flex_vec3_double dpv_dxlo_p(dr_dxlo_p.size());
    for (std::size_t i = 0; i < dr_dxlo_p.size(); ++i) {
      dpv_dxlo_p[i] = crystal_orientation_pv_derivative(
          D, dr_dxlo_p[i], e_X_r, dphi_dxlo_p[i]);
    }
    return dpv_dxlo_p;
  }

  /**
   * Now derivatives of pv and phi wrt each parameter of each crystal unit
   * cell parameterisation that is present.
   */
  vec3 <double> crystal_cell_r_derivative(
      mat3 <double> R, mat3 <double> U, mat3 <double> dB_dxluc_p,
      miller_index h) {
    return R * U * dB_dxluc_p * h;
  }

  /**
   * Now derivatives of r wrt each parameter of each crystal unit
   * cell parameterisation that is present.
   */
  flex_vec3_double crystal_cell_r_derivative(
      mat3 <double> R, mat3 <double> U, const flex_mat3_double &dB_dxluc_p,
      miller_index h) {
    flex_vec3_double dr_dxluc_p(dB_dxluc_p.size());
    for (std::size_t i = 0; i < dB_dxluc_p.size(); ++i) {
      dr_dxluc_p[i] = crystal_cell_r_derivative(R, U, dB_dxluc_p[i], h);
    }
    return dr_dxluc_p;
  }

  /**
   * Now derivatives of phi wrt each parameter of each crystal unit
   * cell parameterisation that is present.
   */
  double crystal_cell_phi_derivative(
      vec3 <double> der, vec3 <double> s, double e_r_s0) {
    return -(der * s) / e_r_s0;
  }

  /**
   * Now derivatives of phi wrt each parameter of each crystal unit
   * cell parameterisation that is present.
   */
  flex_double crystal_cell_phi_derivative(
      const flex_vec3_double &dr_dxluc_p, vec3 <double> s, double e_r_s0) {
    flex_double dphi_dxluc_p(dr_dxluc_p.size());
    for (std::size_t i = 0; i < dr_dxluc_p.size(); ++i) {
      dphi_dxluc_p[i] = crystal_cell_phi_derivative(dr_dxluc_p[i], s, e_r_s0);
    }
    return dphi_dxluc_p;
  }

  /**
   * Now derivatives of pv wrt each parameter of each crystal unit
   * cell parameterisation that is present.
   */
  vec3 <double> crystal_cell_pv_derivative(
      mat3 <double> D, vec3 <double> dr_dxluc_p, vec3 <double> e_X_r,
      double dphi_dxluc_p) {
    return D * (dr_dxluc_p + e_X_r * dphi_dxluc_p);
  }

  /**
   * Now derivatives of pv wrt each parameter of each crystal unit
   * cell parameterisation that is present.
   */
  flex_vec3_double crystal_cell_pv_derivative(
      mat3 <double> D, const flex_vec3_double &dr_dxluc_p, vec3 <double> e_X_r,
      const flex_double &dphi_dxluc_p) {
    SCITBX_ASSERT(dr_dxluc_p.size() == dphi_dxluc_p.size());
    flex_vec3_double dpv_dxluc_p(dr_dxluc_p.size());
    for (std::size_t i = 0; i < dr_dxluc_p.size(); ++i) {
      dpv_dxluc_p[i] = crystal_cell_pv_derivative(
          D, dr_dxluc_p[i], e_X_r, dphi_dxluc_p[i]);
    }
    return dpv_dxluc_p;
  }

}} // namespace dials::refinement

#endif // DIALS_REFINEMENT_PREDICTION_PARAMETER_HELPERS_H
