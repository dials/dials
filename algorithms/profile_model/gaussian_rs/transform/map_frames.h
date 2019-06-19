/*
 * map_frames.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_PROFILE_MODEL_GAUSSIAN_RS_MAP_FRAMES_H
#define DIALS_ALGORITHMS_PROFILE_MODEL_GAUSSIAN_RS_MAP_FRAMES_H

#include <cmath>
#include <algorithm>
#include <boost/math/special_functions/erf.hpp>
#include <scitbx/constants.h>
#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include <scitbx/array_family/tiny_types.h>
#include <dxtbx/model/beam.h>
#include <dxtbx/model/detector.h>
#include <dxtbx/model/goniometer.h>
#include <dxtbx/model/scan.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/error.h>

namespace dials {
  namespace algorithms {
    namespace profile_model {
      namespace gaussian_rs {
  namespace transform {

    using boost::math::erf;
    using scitbx::vec2;

    /**
     * A class to calculate calculate the fraction of counts contributed by each
     * data frame, j, around the reflection to each grid point, v3 in the profile
     * frame.
     */
    template <typename FloatType = double>
    class MapFramesForward {
    public:
      typedef FloatType float_type;

      /**
       * Initialise the fraction calculation class
       * @param starting_angle The angle at frame 0
       * @param oscillation The angular range covered by each frame
       * @param mosaicity The crystal mosaicity
       * @param n_sigma The number of standard deviation to use
       * @param grid_size_e3 The size of the grid
       */
      MapFramesForward(int starting_frame,
                       double starting_angle,
                       double oscillation,
                       double mosaicity,
                       double n_sigma,
                       std::size_t grid_size_e3)
          : starting_frame_(starting_frame),
            starting_angle_(starting_angle),
            oscillation_(oscillation),
            mosaicity_(mosaicity),
            delta_mosaicity_(mosaicity_ * n_sigma),
            grid_size_e3_(grid_size_e3),
            step_size_e3_(2.0 * delta_mosaicity_ / (2 * grid_size_e3_ + 1)) {}

      af::versa<FloatType, af::c_grid<2> > operator()(vec2<int> frames,
                                                      double phi,
                                                      double zeta) const;

    private:
      int starting_frame_;
      double starting_angle_;
      double oscillation_;
      double mosaicity_;
      double delta_mosaicity_;
      std::size_t grid_size_e3_;
      double step_size_e3_;
    };

    /**
     * Calculate the fraction of counts contributed by each data frame, j,
     * around the reflection to each grid point, v3 in the profile frame.
     *
     * First we find and integrate over the set of phi angles covered by each
     * data frame around the reflection to get Ij. Then we find the set of phi
     * angles covered by each grid coordinate. We then integrate over the
     *
     * intersection of the phi angles covered by each data frame and each
     * grid point to get Iv3j. The fraction of the counts is then calculated as
     * Iv3j / Ij.
     *
     * Further details of the method can be found in Kabsch 2010.
     *
     * @param frames The range of frames
     * @param phi The rotation angle of the reflection
     * @param zeta The lorentz correction factor
     *
     * @returns An array containing the count fractions. The fraction of counts
     *          given by frame j to grid coordinate v3 can be found in the array
     *          by fv3j[j-j0, v3-v30]
     *
     * @throws std::runtime_error if the supplied values are bad
     */
    template <typename FloatType>
    af::versa<FloatType, af::c_grid<2> > MapFramesForward<FloatType>::
    operator()(vec2<int> frames, double phi, double zeta) const {
      // Check the value of zeta
      DIALS_ASSERT(frames[1] > frames[0]);
      DIALS_ASSERT(zeta != 0.0);

      // The range of data frames and grid points to iterate over
      int j0 = frames[0];
      int j1 = frames[1];
      int v30 = 0;
      int v31 = 2 * grid_size_e3_ + 1;
      int offset = grid_size_e3_;

      // Create an array to contain the intensity fractions
      af::versa<FloatType, af::c_grid<2> > fraction(
        af::c_grid<2>(j1 - j0, (2 * grid_size_e3_ + 1)),
        af::init_functor_null<FloatType>());

      // A constant used in the solution to the integrals below.
      DIALS_ASSERT(mosaicity_ > 0);
      double sigr2 = std::abs(zeta) / (std::sqrt(2.0) * mosaicity_);

      // Loop over all j data frames in the region around the reflection
      for (int i = 0, j = j0; j < j1; ++j) {
        // The data frame j covers the range of phi such that
        // rj = {phi':phi0 + j*dphi <= phi' >= phi0 + (j+1)*dpi}
        // Therefore the range of phi for j is given as follows.
        double aj = starting_angle_ + (j - starting_frame_) * oscillation_;
        double bj = aj + oscillation_;

        // Calculate the integral over rj (leaving out scaling factors):
        // I[exp(-(phi' - phi)^2 / (2 sigma^2)]
        double integral_j = (erf((bj - phi) * sigr2) - erf((aj - phi) * sigr2));

        // If integral is zero then set fractions to 0.0
        if (integral_j == 0.0) {
          for (int v3 = v30; v3 < v31; ++v3, ++i) {
            fraction[i] = 0.0;
          }
        } else {
          double integral_j_r = 1.0 / integral_j;

          // Loop over all v3 in the profile grid
          for (int v3 = v30; v3 < v31; ++v3) {
            // The grid coordinate v3 cover the range phi such that
            // rv3 = {phi':(v3 - 0.5)d3 <= (phi' - phi)zeta <= (v3 + 0.5)d3}
            // Therefore the range of phi for v3 is given as follows.
            double bv3 = ((v3 - offset - 0.5) * step_size_e3_) / zeta + phi;
            double av3 = ((v3 - offset + 0.5) * step_size_e3_) / zeta + phi;
            if (av3 > bv3) std::swap(av3, bv3);

            // We need to integrate over the intersection of sets rv3 and rj
            double av3j = std::max(av3, aj);
            double bv3j = std::min(bv3, bj);

            // If there is no intersection then set the fraction of the
            // counts contributed by data frame j to grid coordinate v3 to
            // zero, otherwise calculate it as the ratio of the integral
            // over the intersection or rv3 and rj to the integral over rj
            if (av3j >= bv3j) {
              fraction[i] = 0.0;
            } else {
              fraction[i] = (FloatType)(
                (erf((bv3j - phi) * sigr2) - erf((av3j - phi) * sigr2)) * integral_j_r);
            }

            // Increment array index
            i++;
          }
        }
      }

      // Return the intensity fractions
      return fraction;
    }

    /**
     * A class to calculate calculate the fraction of counts contributed by each
     * grid point, v3 to each data frame, j in the profile.
     */
    template <typename FloatType = double>
    class MapFramesReverse {
    public:
      typedef FloatType float_type;

      /**
       * Initialise the fraction calculation class
       * @param starting_angle The angle at frame 0
       * @param oscillation The angular range covered by each frame
       * @param mosaicity The crystal mosaicity
       * @param n_sigma The number of standard deviation to use
       * @param grid_size_e3 The size of the grid
       */
      MapFramesReverse(int starting_frame,
                       double starting_angle,
                       double oscillation,
                       double mosaicity,
                       double n_sigma,
                       std::size_t grid_size_e3)
          : starting_frame_(starting_frame),
            starting_angle_(starting_angle),
            oscillation_(oscillation),
            mosaicity_(mosaicity),
            delta_mosaicity_(mosaicity_ * n_sigma),
            grid_size_e3_(grid_size_e3),
            step_size_e3_(2.0 * delta_mosaicity_ / (2 * grid_size_e3_ + 1)) {}

      af::versa<FloatType, af::c_grid<2> > operator()(vec2<int> bbox_z,
                                                      double phi,
                                                      double zeta) const;

    private:
      int starting_frame_;
      double starting_angle_;
      double oscillation_;
      double mosaicity_;
      double delta_mosaicity_;
      std::size_t grid_size_e3_;
      double step_size_e3_;
    };

    /**
     * Calculate the fraction of counts contributed by each grid point, v3,
     * to each data frame, j.
     *
     * First we find and integrate over the set of phi angles covered by each
     * grid point in the profile to get Iv3. Then we find the set of phi
     * angles covered by each data frame. We then integrate over the
     * intersection of the phi angles covered by each data frame and each
     * grid point to get Iv3j. The fraction of the counts is then calculated as
     * Iv3j / Iv3.
     *
     * This is the reverse of the forward procedure.
     *
     * @param frames The range of frames
     * @param phi The rotation angle of the reflection
     * @param zeta The lorentz correction factor
     *
     * @returns An array containing the count fractions. The fraction of counts
     *          given by grid point v3 to data frame j can be found in the array
     *          by fv3j[v3-v30, j-j0]
     *
     * @throws std::runtime_error if the supplied values are bad
     */
    template <typename FloatType>
    af::versa<FloatType, af::c_grid<2> > MapFramesReverse<FloatType>::
    operator()(vec2<int> frames, double phi, double zeta) const {
      // Check the value of zeta
      DIALS_ASSERT(frames[1] > frames[0]);
      DIALS_ASSERT(zeta != 0.0);

      // The range of data frames and grid points to iterate over
      int j0 = frames[0];
      int j1 = frames[1];
      int v30 = -(int)grid_size_e3_;
      int v31 = +(int)grid_size_e3_ + 1;

      // Create an array to contain the intensity fractions
      af::versa<FloatType, af::c_grid<2> > fraction(
        af::c_grid<2>((2 * grid_size_e3_ + 1), j1 - j0),
        af::init_functor_null<FloatType>());

      // A constant used in the solution to the integrals below.
      DIALS_ASSERT(mosaicity_ > 0);
      double sigr2 = std::abs(zeta) / (std::sqrt(2.0) * mosaicity_);

      // Loop over all v3 grid points in the profile
      for (int i = 0, v3 = v30; v3 < v31; ++v3) {
        // The grid coordinate v3 cover the range phi such that
        // rv3 = {phi':(v3 - 0.5)d3 <= (phi' - phi)zeta <= (v3 + 0.5)d3}
        // Therefore the range of phi for v3 is given as follows.
        double bv3 = ((v3 - 0.5) * step_size_e3_) / zeta + phi;
        double av3 = ((v3 + 0.5) * step_size_e3_) / zeta + phi;
        if (av3 > bv3) std::swap(av3, bv3);

        // Calculate the integral over rv3 (leaving out scaling factors):
        // I[exp(-(phi' - phi)^2 / (2 sigma^2)]
        double integral_v3 = (erf((bv3 - phi) * sigr2) - erf((av3 - phi) * sigr2));

        // If integral is zero then set fractions to 0.0
        if (integral_v3 == 0.0) {
          for (int j = j0; j < j1; ++j, ++i) {
            fraction[i] = 0.0;
          }
        } else {
          double integral_v3_r = 1.0 / integral_v3;

          // Loop over all data frames
          for (int j = j0; j < j1; ++j) {
            // The data frame j covers the range of phi such that
            // rj = {phi':phi0 + j*dphi <= phi' >= phi0 + (j+1)*dpi}
            // Therefore the range of phi for j is given as follows.
            double aj = starting_angle_ + (j - starting_frame_) * oscillation_;
            double bj = aj + oscillation_;

            // We need to integrate over the intersection of sets rv3 and rj
            double av3j = std::max(av3, aj);
            double bv3j = std::min(bv3, bj);

            // If there is no intersection then set the fraction of the
            // counts contributed by data frame j to grid coordinate v3 to
            // zero, otherwise calculate it as the ratio of the integral
            // over the intersection or rv3 and rj to the integral over rj
            if (av3j >= bv3j) {
              fraction[i] = 0.0;
            } else {
              fraction[i] =
                (FloatType)((erf((bv3j - phi) * sigr2) - erf((av3j - phi) * sigr2))
                            * integral_v3_r);
            }

            // Increment array index
            i++;
          }
        }
      }

      // Return the intensity fractions
      return fraction;
    }

}}}}}  // namespace dials::algorithms::profile_model::gaussian_rs::transform

#endif  // DIALS_ALGORITHMS_PROFILE_MODEL_GAUSSIAN_RS_MAP_FRAMES_H
