/*
 * corrections.cc
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */

#ifndef DIALS_ALGORITHMS_INTEGRATION_CORRECTIONS_H
#define DIALS_ALGORITHMS_INTEGRATION_CORRECTIONS_H

#include <vector>
#include <scitbx/vec3.h>
#include <dxtbx/model/beam.h>
#include <dxtbx/model/detector.h>
#include <dxtbx/model/goniometer.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  using dxtbx::model::BeamBase;
  using dxtbx::model::Detector;
  using dxtbx::model::Goniometer;
  using scitbx::vec3;

  /**
   * Compute the LP correction for a single reflection. Note that the
   * polarization fraction follows XDS convention, in which a value of 0.5
   * implies an unpolarized beam, rather than the MOSFLM definition in which
   * an unpolarized beam has a polarization factor of 0.0.
   * @param s0 The incident beam vector
   * @param pn The polarization plane normal
   * @param pf The polarization plane fraction
   * @param m2 The rotation axis
   * @param s1 The incident beam vector
   * @returns L / P The correction factor
   */
  double lp_correction(vec3<double> s0,
                       vec3<double> pn,
                       double pf,
                       vec3<double> m2,
                       vec3<double> s1) {
    double s1_length = s1.length();
    double s0_length = s0.length();
    DIALS_ASSERT(s1_length > 0 && s0_length > 0);
    double L = std::abs(s1 * (m2.cross(s0))) / (s0_length * s1_length);
    double P1 = ((pn * s1) / s1_length);
    double P2 = (1.0 - 2.0 * pf) * (1.0 - P1 * P1);
    double P3 = (s1 * s0 / (s1_length * s0_length));
    double P4 = pf * (1.0 + P3 * P3);
    double P = P2 + P4;
    DIALS_ASSERT(P != 0);
    return L / P;
  }

  /**
   * Compute the Stills LP correction for a single reflection. Note that the
   * polarization fraction follows XDS convention, in which a value of 0.5
   * implies an unpolarized beam, rather than the MOSFLM definition in which
   * an unpolarized beam has a polarization factor of 0.0.
   * @param s0 The incident beam vector
   * @param pn The polarization plane normal
   * @param pf The polarization plane fraction
   * @param s1 The incident beam vector
   * @returns L / P The correction factor
   */
  double stills_lp_correction(vec3<double> s0,
                              vec3<double> pn,
                              double pf,
                              vec3<double> s1) {
    double s1_length = s1.length();
    double s0_length = s0.length();
    DIALS_ASSERT(s1_length > 0 && s0_length > 0);
    double P1 = ((pn * s1) / s1_length);
    double P2 = (1.0 - 2.0 * pf) * (1.0 - P1 * P1);
    double P3 = (s1 * s0 / (s1_length * s0_length));
    double P4 = pf * (1.0 + P3 * P3);
    double P = P2 + P4;
    DIALS_ASSERT(P != 0);
    return 1.0 / P;
  }

  /**
   * Compute the QE correction for a single reflection
   * @param mu attenuation coefficient in mm^-1
   * @param t0 thickness of sensor in mm
   * @param s1 direction of diffracted ray
   * @param n detector / panel normal for this reflection
   * @returns QE term which needs to be divided by (i.e. is efficiency)
   */

  double qe_correction(double mu, double t0, vec3<double> s1, vec3<double> n) {
    DIALS_ASSERT(mu >= 0);
    DIALS_ASSERT(t0 >= 0);
    double cos_angle = cos(n.angle(s1));
    cos_angle = std::abs(cos_angle);
    double t = t0 / cos_angle;
    return 1.0 - exp(-mu * t);
  }

  /**
   * A class to perform corrections to the intensities.
   */
  class Corrections {
  public:
    /**
     * @param beam The beam model.
     * @param goniometer The goniometer model.
     */
    Corrections(const BeamBase &beam, const Goniometer &goniometer)
        : s0_(beam.get_s0()),
          pn_(beam.get_polarization_normal()),
          pf_(beam.get_polarization_fraction()),
          m2_(goniometer.get_rotation_axis()) {
      // Deprecated constructor
    }

    /**
     * @param beam The beam model.
     * @param goniometer The goniometer model.
     * @param detector The detector model.
     */
    Corrections(const BeamBase &beam,
                const Goniometer &goniometer,
                const Detector &detector)
        : s0_(beam.get_s0()),
          pn_(beam.get_polarization_normal()),
          pf_(beam.get_polarization_fraction()),
          m2_(goniometer.get_rotation_axis()),
          det_(detector) {}

    /**
     * @param beam The beam model.
     * @param goniometer The goniometer model.
     * @param detector The detector model.
     */
    Corrections(const BeamBase &beam, const Detector &detector)
        : s0_(beam.get_s0()),
          pn_(beam.get_polarization_normal()),
          pf_(beam.get_polarization_fraction()),
          m2_(0, 0, 0),
          det_(detector) {}

    /**
     * Perform the LP correction. If no rotation axis is specified then do the
     * stills lorentz correction
     * @param s1 The incident beam vector
     * @returns L / P The correction
     */
    double lp(vec3<double> s1) const {
      if (m2_.length() > 0) {
        return lp_correction(s0_, pn_, pf_, m2_, s1);
      } else {
        return stills_lp_correction(s0_, pn_, pf_, s1);
      }
    }

    /**
     * Perform the QE correction
     * @param s1 The incident beam vector
     * @param p The panel for this reflection
     * @returns QE term which needs to be divided by (i.e. is efficiency)
     */
    double qe(vec3<double> s1, size_t p) const {
      return qe_correction(
        det_[p].get_mu(), det_[p].get_thickness(), s1, det_[p].get_normal());
    }

  private:
    vec3<double> s0_;
    vec3<double> pn_;
    double pf_;
    vec3<double> m2_;
    Detector det_;
  };

  /**
   * A class to perform corrections for multiple experiments
   */
  class CorrectionsMulti {
  public:
    /**
     * Add another correction class
     * @param obj The correction class
     */
    void push_back(const Corrections &obj) {
      compute_.push_back(obj);
    }

    /**
     * @returns The number of correctors.
     */
    std::size_t size() const {
      return compute_.size();
    }

    /**
     * Perform the LP correction.
     * @param id The list of experiments ids
     * @param s1 The list of incident beam vectors
     */
    af::shared<double> lp(const af::const_ref<int> &id,
                          const af::const_ref<vec3<double> > &s1) const {
      DIALS_ASSERT(id.size() == s1.size());
      af::shared<double> result(id.size(), 0);
      for (std::size_t i = 0; i < id.size(); ++i) {
        DIALS_ASSERT(id[i] >= 0);
        DIALS_ASSERT(id[i] < compute_.size());
        result[i] = compute_[id[i]].lp(s1[i]);
      }
      return result;
    }

    /**
     * Perform the QE correction.
     * @param id The list of experiments ids
     * @param s1 The list of incident beam vectors
     * @param p The list of panels
     */
    af::shared<double> qe(const af::const_ref<int> &id,
                          const af::const_ref<vec3<double> > &s1,
                          const af::const_ref<std::size_t> &p) const {
      DIALS_ASSERT(id.size() == s1.size());
      DIALS_ASSERT(id.size() == p.size());
      af::shared<double> result(id.size(), 0);
      for (std::size_t i = 0; i < id.size(); ++i) {
        DIALS_ASSERT(id[i] >= 0);
        DIALS_ASSERT(id[i] < compute_.size());
        result[i] = compute_[id[i]].qe(s1[i], p[i]);
      }
      return result;
    }

  private:
    std::vector<Corrections> compute_;
  };
}}  // namespace dials::algorithms

#endif  // DIALS_ALGORITHMS_INTEGRATION_CORRECTIONS_H
