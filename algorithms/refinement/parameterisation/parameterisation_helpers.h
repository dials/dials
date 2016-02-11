
#ifndef DIALS_REFINEMENT_PREDICTION_PARAMETER_HELPERS_H
#define DIALS_REFINEMENT_PREDICTION_PARAMETER_HELPERS_H

#include <scitbx/vec3.h>
#include <scitbx/mat3.h>
#include <cctbx/miller.h>
#include <dxtbx/model/detector.h>
#include <dxtbx/model/panel.h>
#include <dials/error.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
//#include <dials/algorithms/refinement/rtmats.h>
#include <scitbx/math/r3_rotation.h>

namespace dials { namespace refinement {

  using scitbx::vec3;
  using scitbx::mat3;
  using dxtbx::model::Detector;
  using scitbx::math::r3_rotation::axis_and_angle_as_matrix;

  // declare this here - can't #include the header as the function is defined
  // there, leading to multiple definitions
  mat3<double> dR_from_axis_and_angle(const vec3<double> &axis,
                                      double angle,
                                      bool deg=false);

  af::shared < mat3 <double> >
  selected_multi_panel_compose(
      const af::const_ref< vec3 <double> > &initial_state,
      const af::const_ref< double> &params_vals,
      const af::const_ref< vec3 <double> > &params_axes,
      Detector &detector,
      const af::const_ref<int> &selection,
      const af::const_ref< vec3 <double> > &offsets,
      const af::const_ref< vec3 <double> > &dir1s,
      const af::const_ref< vec3 <double> > &dir2s,
      const mat3 <double> & Tau1, const mat3 <double> & dTau1_dtau1,
      const mat3 <double> & Tau2, const mat3 <double> & dTau2_dtau2,
      const mat3 <double> & Tau3, const mat3 <double> & dTau3_dtau3) {

    af::shared < mat3 <double> >  ret(6*selection.size(),af::init_functor_null< mat3<double> >());

    // extract items from the initial state
    vec3<double> id1 = initial_state[0];
    vec3<double> id2 = initial_state[1];

    // extract parameters from the internal list.  v: value, a: axis
    double dist_v   = params_vals[0];
    double shift1_v = params_vals[1];
    double shift2_v = params_vals[2];

    vec3<double> dist_a   = params_axes[0];
    vec3<double> shift1_a = params_axes[1];
    vec3<double> shift2_a = params_axes[2];

    // finish composing derivatives
    mat3 <double> Tau32 = Tau3 * Tau2;
    mat3 <double> Tau321 = Tau32 * Tau1;

    /*
    # Compose new state
    # =================
    # First the frame positioned at a distance from the lab origin
    */
    vec3<double> P0 = dist_v * dist_a; // distance along initial detector normal
    vec3<double> Px = P0 + id1;        // point at the end of d1 in lab frame
    vec3<double> Py = P0 + id2;        // point at the end of d2 in lab frame

    // detector shift vector
    vec3<double> dsv = P0 + shift1_v * shift1_a + shift2_v * shift2_a;

    // compose dorg point
    vec3<double> dorg = Tau321 * dsv - Tau32 * P0 + P0;

    // compose new d1, d2 and dn and ensure frame remains orthonormal.
    vec3<double> d1 = (Tau321 * (Px - P0)).normalize();
    vec3<double> d2 = (Tau321 * (Py - P0)).normalize();
    vec3<double> dn = d1.cross(d2).normalize();
    d2 = dn.cross(d1);

    // compose new Panel origins
    af::shared< vec3<double> > origins(offsets.size(),
      af::init_functor_null< vec3<double> >());
    for(std::size_t i = 0; i < offsets.size(); i++)
      origins[i] = dorg + offsets[i][0] * d1 + \
                          offsets[i][1] * d2 + \
                          offsets[i][2] * dn;

    // compose new Panel directions
    af::shared< vec3<double> > dir1s_new(dir1s.size(),
      af::init_functor_null< vec3<double> >());
    for(std::size_t i = 0; i < dir1s.size(); i++)
      dir1s_new[i] = dir1s[i][0] * d1 + \
                     dir1s[i][1] * d2 + \
                     dir1s[i][2] * dn;

    af::shared< vec3<double> > dir2s_new(dir2s.size(),
      af::init_functor_null< vec3<double> >());
    for(std::size_t i = 0; i < dir2s.size(); i++)
      dir2s_new[i] = dir2s[i][0] * d1 + \
                     dir2s[i][1] * d2 + \
                     dir2s[i][2] * dn;

    DIALS_ASSERT(selection.size() == dir1s_new.size() &&
                 selection.size() == dir2s_new.size() &&
                 selection.size() == origins.size());

    // now update the panels with their new position and orientation.
    for(std::size_t i = 0; i < selection.size(); i++)
      detector[selection[i]].set_frame(dir1s_new[i], dir2s_new[i], origins[i]);
    /*
    # calculate derivatives of the state wrt parameters
    # =================================================
    # Start with the dorg vector, where
    # dorg = Tau321 * dsv - Tau32 * P0 + P0
    */

    // derivative wrt dist
    vec3 <double> dP0_ddist   = dist_a;
    vec3 <double> ddsv_ddist  = dP0_ddist;
    vec3 <double> ddorg_ddist = Tau321 * ddsv_ddist - Tau32 * dP0_ddist + dP0_ddist;

    // derivative wrt shift1
    vec3 <double> ddsv_dshift1  = shift1_a;
    vec3 <double> ddorg_dshift1 = Tau321 * ddsv_dshift1;

    // derivative wrt shift2
    vec3 <double> ddsv_dshift2  = shift2_a;
    vec3 <double> ddorg_dshift2 = Tau321 * ddsv_dshift2;

    // derivative wrt tau1
    mat3 <double> dTau321_dtau1 = Tau32 * dTau1_dtau1;
    vec3 <double> ddorg_dtau1   = dTau321_dtau1 * dsv;

    // derivative wrt tau2
    mat3 <double> dTau32_dtau2  = Tau3 * dTau2_dtau2;
    mat3 <double> dTau321_dtau2 = dTau32_dtau2 * Tau1;
    vec3 <double> ddorg_dtau2   = dTau321_dtau2 * dsv - dTau32_dtau2 * P0;

    // derivative wrt tau3
    mat3 <double> dTau32_dtau3 = dTau3_dtau3 * Tau2;
    mat3 <double> dTau321_dtau3 = dTau32_dtau3 * Tau1;
    vec3 <double> ddorg_dtau3 = dTau321_dtau3 * dsv - dTau32_dtau3 * P0;

    /*
    # Now derivatives of the direction d1, where
    # d1 = (Tau321 * (Px - P0)).normalize()
    # For calc of derivatives ignore the normalize(), which should
    # be unnecessary anyway as Px - P0 is a unit vector and Tau321 a
    # pure rotation.
    */

    // derivative wrt dist
    // dPx_ddist = dist.axis; dP0_ddist = dist.axis, so these cancel
    vec3 <double> dd1_ddist(0., 0., 0.);

    // derivative wrt shift1
    vec3 <double> dd1_dshift1(0., 0., 0.);

    // derivative wrt shift2
    vec3 <double> dd1_dshift2(0., 0., 0.);

    // derivative wrt tau1
    vec3 <double> dd1_dtau1 = dTau321_dtau1 * (Px - P0);

    // derivative wrt tau2
    vec3 <double> dd1_dtau2 = dTau321_dtau2 * (Px - P0);

    // derivative wrt tau3
    vec3 <double> dd1_dtau3 = dTau321_dtau3 * (Px - P0);

    // Derivatives of the direction d2, where
    // d2 = (Tau321 * (Py - P0)).normalize()

    // derivative wrt dist
    vec3 <double> dd2_ddist(0., 0., 0.);

    // derivative wrt shift1
    vec3 <double> dd2_dshift1(0., 0., 0.);

    // derivative wrt shift2
    vec3 <double> dd2_dshift2(0., 0., 0.);

    // derivative wrt tau1
    vec3 <double> dd2_dtau1 = dTau321_dtau1 * (Py - P0);

    // derivative wrt tau2
    vec3 <double> dd2_dtau2 = dTau321_dtau2 * (Py - P0);

    // derivative wrt tau3
    vec3 <double> dd2_dtau3 = dTau321_dtau3 * (Py - P0);

    // Derivatives of the direction dn, where
    // dn = d1.cross(d2).normalize()

    // derivative wrt dist
    vec3 <double> ddn_ddist(0., 0., 0.);

    // derivative wrt shift1
    vec3 <double> ddn_dshift1(0., 0., 0.);

    // derivative wrt shift2
    vec3 <double> ddn_dshift2(0., 0., 0.);

    // derivative wrt tau1. Product rule for cross product applies
    vec3 <double> ddn_dtau1 = dd1_dtau1.cross(d2) + d1.cross(dd2_dtau1);

    // derivative wrt tau2
    vec3 <double> ddn_dtau2 = dd1_dtau2.cross(d2) + d1.cross(dd2_dtau2);

    // derivative wrt tau3
    vec3 <double> ddn_dtau3 = dd1_dtau3.cross(d2) + d1.cross(dd2_dtau3);

    // calculate derivatives of the attached Panel matrices
    //====================================================
    for (std::size_t sel_id = 0; sel_id < selection.size(); sel_id++) {
      vec3 <double> offset = offsets[sel_id];
      vec3 <double> dir1_new_basis = dir1s[sel_id];
      vec3 <double> dir2_new_basis = dir2s[sel_id];

      // Panel origin:
      // o = dorg + offset[0] * d1 + offset[1] * d2 + offset[2] * dn

      // derivative wrt dist. NB only ddorg_ddist is not null! The other
      // elements are left here to aid understanding, but should be removed
      // when this class is ported to C++ for speed.
      vec3 <double> do_ddist = ddorg_ddist + offset[0] * dd1_ddist + \
                                             offset[1] * dd2_ddist + \
                                             offset[2] * ddn_ddist;

      // derivative wrt shift1. NB only ddorg_dshift1 is non-null.
      vec3 <double> do_dshift1 = ddorg_dshift1 + offset[0] * dd1_dshift1 + \
                                                 offset[1] * dd2_dshift1 + \
                                                 offset[2] * ddn_dshift1;

      // derivative wrt shift2. NB only ddorg_dshift2 is non-null.
      vec3 <double> do_dshift2 = ddorg_dshift2 + offset[0] * dd1_dshift2 + \
                                                 offset[1] * dd2_dshift2 + \
                                                 offset[2] * ddn_dshift2;

      // derivative wrt tau1
      vec3 <double> do_dtau1 = ddorg_dtau1 + offset[0] * dd1_dtau1 + \
                                             offset[1] * dd2_dtau1 + \
                                             offset[2] * ddn_dtau1;

      // derivative wrt tau2
      vec3 <double> do_dtau2 = ddorg_dtau2 + offset[0] * dd1_dtau2 + \
                                             offset[1] * dd2_dtau2 + \
                                             offset[2] * ddn_dtau2;

      // derivative wrt tau3
      vec3 <double> do_dtau3 = ddorg_dtau3 + offset[0] * dd1_dtau3 + \
                                             offset[1] * dd2_dtau3 + \
                                             offset[2] * ddn_dtau3;

      // Panel dir1:
      // dir1 = dir1_new_basis[0] * d1 + dir1_new_basis[1] * d2 +
      //        dir1_new_basis[2] * dn

      // derivative wrt dist. NB These are all null.
      vec3 <double> ddir1_ddist = dir1_new_basis[0] * dd1_ddist + \
                                  dir1_new_basis[1] * dd2_ddist + \
                                  dir1_new_basis[2] * ddn_ddist;

      // derivative wrt shift1. NB These are all null.
      vec3 <double> ddir1_dshift1 = dir1_new_basis[0] * dd1_dshift1 + \
                                    dir1_new_basis[1] * dd2_dshift1 + \
                                    dir1_new_basis[2] * ddn_dshift1;

      // derivative wrt shift2. NB These are all null.
      vec3 <double> ddir1_dshift2 = dir1_new_basis[0] * dd1_dshift2 + \
                                    dir1_new_basis[1] * dd2_dshift2 + \
                                    dir1_new_basis[2] * ddn_dshift2;

      // derivative wrt tau1
      vec3 <double> ddir1_dtau1 = dir1_new_basis[0] * dd1_dtau1 + \
                                  dir1_new_basis[1] * dd2_dtau1 + \
                                  dir1_new_basis[2] * ddn_dtau1;

      // derivative wrt tau2
      vec3 <double> ddir1_dtau2 = dir1_new_basis[0] * dd1_dtau2 + \
                                  dir1_new_basis[1] * dd2_dtau2 + \
                                  dir1_new_basis[2] * ddn_dtau2;

      // derivative wrt tau3
      vec3 <double> ddir1_dtau3 = dir1_new_basis[0] * dd1_dtau3 + \
                                  dir1_new_basis[1] * dd2_dtau3 + \
                                  dir1_new_basis[2] * ddn_dtau3;

      // Panel dir2:
      // dir2 = dir2_new_basis[0] * d1 + dir2_new_basis[1] * d2 +
      //        dir2_new_basis[2] * dn

      // derivative wrt dist. NB These are all null.
      vec3 <double> ddir2_ddist = dir2_new_basis[0] * dd1_ddist + \
                                  dir2_new_basis[1] * dd2_ddist + \
                                  dir2_new_basis[2] * ddn_ddist;

      // derivative wrt shift1. NB These are all null.
      vec3 <double> ddir2_dshift1 = dir2_new_basis[0] * dd1_dshift1 + \
                                    dir2_new_basis[1] * dd2_dshift1 + \
                                    dir2_new_basis[2] * ddn_dshift1;

      // derivative wrt shift2. NB These are all null.
      vec3 <double> ddir2_dshift2 = dir2_new_basis[0] * dd1_dshift2 + \
                                    dir2_new_basis[1] * dd2_dshift2 + \
                                    dir2_new_basis[2] * ddn_dshift2;

      // derivative wrt tau1
      vec3 <double> ddir2_dtau1 = dir2_new_basis[0] * dd1_dtau1 + \
                                  dir2_new_basis[1] * dd2_dtau1 + \
                                  dir2_new_basis[2] * ddn_dtau1;

      // derivative wrt tau2
      vec3 <double> ddir2_dtau2 = dir2_new_basis[0] * dd1_dtau2 + \
                                  dir2_new_basis[1] * dd2_dtau2 + \
                                  dir2_new_basis[2] * ddn_dtau2;

      // derivative wrt tau3
      vec3 <double> ddir2_dtau3 = dir2_new_basis[0] * dd1_dtau3 + \
                                  dir2_new_basis[1] * dd2_dtau3 + \
                                  dir2_new_basis[2] * ddn_dtau3;

      // combine these vectors together into derivatives of the panel
      // matrix d and return them, converting angles back to mrad

      // derivative wrt dist
      ret[0*selection.size() + sel_id] = mat3<double>(ddir1_ddist[0],ddir1_ddist[1],ddir1_ddist[2],
                                                      ddir2_ddist[0],ddir2_ddist[1],ddir2_ddist[2],
                                                      do_ddist[0],   do_ddist[1],   do_ddist[2]).transpose();

      // derivative wrt shift1
      ret[1*selection.size() + sel_id] = mat3<double>(ddir1_dshift1[0],ddir1_dshift1[1],ddir1_dshift1[2],
                                                      ddir2_dshift1[0],ddir2_dshift1[1],ddir2_dshift1[2],
                                                      do_dshift1[0],   do_dshift1[1],   do_dshift1[2]).transpose();

      // derivative wrt shift2
      ret[2*selection.size() + sel_id] = mat3<double>(ddir1_dshift2[0],ddir1_dshift2[1],ddir1_dshift2[2],
                                                      ddir2_dshift2[0],ddir2_dshift2[1],ddir2_dshift2[2],
                                                      do_dshift2[0],   do_dshift2[1],   do_dshift2[2]).transpose();

      // derivative wrt tau1
      ret[3*selection.size() + sel_id] = mat3<double>(ddir1_dtau1[0],ddir1_dtau1[1],ddir1_dtau1[2],
                                                      ddir2_dtau1[0],ddir2_dtau1[1],ddir2_dtau1[2],
                                                      do_dtau1[0],   do_dtau1[1],   do_dtau1[2]).transpose() / 1000.;

      // derivative wrt tau2
      ret[4*selection.size() + sel_id] = mat3<double>(ddir1_dtau2[0],ddir1_dtau2[1],ddir1_dtau2[2],
                                                      ddir2_dtau2[0],ddir2_dtau2[1],ddir2_dtau2[2],
                                                      do_dtau2[0],   do_dtau2[1],   do_dtau2[2]).transpose() / 1000.;

      // derivative wrt tau3
      ret[5*selection.size() + sel_id] = mat3<double>(ddir1_dtau3[0],ddir1_dtau3[1],ddir1_dtau3[2],
                                                      ddir2_dtau3[0],ddir2_dtau3[1],ddir2_dtau3[2],
                                                      do_dtau3[0],   do_dtau3[1],   do_dtau3[2]).transpose() / 1000.;
    }
    return ret;
  }

  af::shared < mat3 <double> >
  multi_panel_compose(
      const af::const_ref< vec3 <double> > &initial_state,
      const af::const_ref< double> &params_vals,
      const af::const_ref< vec3 <double> > &params_axes,
      Detector &detector,
      const af::const_ref< vec3 <double> > &offsets,
      const af::const_ref< vec3 <double> > &dir1s,
      const af::const_ref< vec3 <double> > &dir2s,
      const mat3 <double> & Tau1, const mat3 <double> & dTau1_dtau1,
      const mat3 <double> & Tau2, const mat3 <double> & dTau2_dtau2,
      const mat3 <double> & Tau3, const mat3 <double> & dTau3_dtau3) {

    af::shared<int> selection;
    selection.reserve(detector.size());
    for (int i = 0; i < detector.size(); i++) {
      selection.push_back(i);
    }
    return selected_multi_panel_compose(initial_state, params_vals, params_axes,
      detector, selection.const_ref(), offsets, dir1s, dir2s, Tau1, dTau1_dtau1,
      Tau2, dTau2_dtau2, Tau3, dTau3_dtau3);
  }

  /**
   * Given an initial orientation matrix, the values of the three orientation
   * parameters, phi1, phi2 and phi3 (in mrad), plus the axes about which these
   * parameters act, calculate the new orientation U plus the derivatives
   * dU/dphi1, dU/dphi2 and dU/dphi3
   */
  class CrystalOrientationCompose {

  public:

    CrystalOrientationCompose(
        const mat3<double> &U0,
        double phi1,
        const vec3<double> &phi1_axis,
        double phi2,
        const vec3<double> &phi2_axis,
        double phi3,
        const vec3<double> &phi3_axis) {

      // convert angles from mrad to radians
      phi1 /= 1000.;
      phi2 /= 1000.;
      phi3 /= 1000.;

      // compose rotation matrices and their first order derivatives
      mat3<double> Phi1 = axis_and_angle_as_matrix(phi1_axis, phi1);
      mat3<double> dPhi1_dphi1 = dR_from_axis_and_angle(phi1_axis, phi1);

      mat3<double> Phi2 = axis_and_angle_as_matrix(phi2_axis, phi2);
      mat3<double> dPhi2_dphi2 = dR_from_axis_and_angle(phi2_axis, phi2);

      mat3<double> Phi3 = axis_and_angle_as_matrix(phi3_axis, phi3);
      mat3<double> dPhi3_dphi3 = dR_from_axis_and_angle(phi3_axis, phi3);

      // compose new state
      mat3<double> Phi21 = Phi2 * Phi1;
      mat3<double> Phi321 = Phi3 * Phi21;
      U_ = Phi321 * U0;

      // calculate derivatives of the state wrt parameters
      dU_dphi1_ = (Phi3 * Phi2 * dPhi1_dphi1 * U0)/1000.0;
      dU_dphi2_ = (Phi3 * dPhi2_dphi2 * Phi1 * U0)/1000.0;
      dU_dphi3_ = (dPhi3_dphi3 * Phi21 * U0)/1000.0;

    }

    mat3<double>
    U(){
      return U_;
    }

    mat3<double>
    dU_dphi1(){
      return dU_dphi1_;
    }

    mat3<double>
    dU_dphi2(){
      return dU_dphi2_;
    }

    mat3<double>
    dU_dphi3(){
      return dU_dphi3_;
    }

  private:
    mat3<double> U_;
    mat3<double> dU_dphi1_;
    mat3<double> dU_dphi2_;
    mat3<double> dU_dphi3_;
  };

}} // namespace dials::refinement

#endif // DIALS_REFINEMENT_PREDICTION_PARAMETER_HELPERS_H
