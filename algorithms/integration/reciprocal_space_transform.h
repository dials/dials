/*
 * reciprocal_space_transform.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_INTEGRATION_XDS_TRANSFORM_H
#define DIALS_ALGORITHMS_INTEGRATION_XDS_TRANSFORM_H

#include <cmath>
#include <algorithm>
#include <boost/math/special_functions/erf.hpp>
#include <scitbx/constants.h>
#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/tiny_types.h>
#include <dxtbx/model/beam.h>
#include <dxtbx/model/detector.h>
#include <dxtbx/model/goniometer.h>
#include <dxtbx/model/scan.h>
#include <dials/model/data/reflection.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  using boost::math::erf;
  using scitbx::vec2;
  using scitbx::vec3;
  using scitbx::af::int6;
  using scitbx::af::flex_int;
  using scitbx::af::flex_double;
  using scitbx::af::flex_grid;
  using dxtbx::model::Beam;
  using dxtbx::model::Detector;
  using dxtbx::model::Goniometer;
  using dxtbx::model::Scan;
  using dials::model::Reflection;
  using dials::model::ReflectionList;

  typedef scitbx::af::flex<vec2<double> >::type flex_vec2_double;
  typedef scitbx::af::flex<vec3<double> >::type flex_vec3_double;

  /**
   * A class to calculate calculate the fraction of counts contributed by each
   * data frame, j, around the reflection to each grid point, v3 in the profile
   * frame.
   */
  class ReciprocalSpaceTransformE3Fraction {

  public:

    /**
     * Initialise the fraction calculation class
     * @param scan The scan parameters
     * @param mosaicity The crystal mosaicity
     * @param n_sigma The number of standard deviation to use
     * @param grid_size_e3 The size of the grid
     */
    ReciprocalSpaceTransformE3Fraction(const Scan &scan,
                                       double mosaicity,
                                       double n_sigma,
                                       std::size_t grid_size_e3)
      : starting_angle_(scan.get_oscillation()[0]),
        oscillation_(scan.get_oscillation()[1]),
        mosaicity_(mosaicity),
        delta_mosaicity_(mosaicity_ * n_sigma),
        grid_size_e3_(grid_size_e3),
        step_size_e3_(delta_mosaicity_ / (2 * grid_size_e3_ + 1)) {}

    flex_double operator()(vec2 <int> bbox_z, double phi, double zeta) const;

  private:

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
   * @param bbox_z The z bounds of the reflection
   * @param phi The rotation angle of the reflection
   * @param zeta The lorentz correction factor
   *
   * @returns An array containing the count fractions. The fraction of counts
   *          given by frame j to grid coordinate v3 can be found in the array
   *          by fv3j[v3-v30, j-j0]
   *
   * @throws std::runtime_error if the supplied values are bad
   */
  flex_double ReciprocalSpaceTransformE3Fraction::operator()(
      vec2 <int> bbox_z, double phi, double zeta) const
  {
    // Check the value of zeta
    DIALS_ASSERT(bbox_z[0] >= 0 && bbox_z[1] > bbox_z[0]);
    DIALS_ASSERT(zeta != 0.0);

    // The range of data frames and grid points to iterate over
    std::size_t j0 = bbox_z[0];
    std::size_t j1 = bbox_z[1];
    int v30 = - (int)grid_size_e3_;
    int v31 = + (int)grid_size_e3_ + 1;

    // Create an array to contain the intensity fractions
    flex_double fraction(flex_grid<>((2 * grid_size_e3_ + 1), j1 - j0));

    // A constant used in the solution to the integrals below.
    double sigr2 = 1.0 / (std::sqrt(2.0) * (mosaicity_ / std::abs(zeta)));

    // Loop over all j data frames in the region around the reflection
    for (std::size_t i = 0, j = j0; j < j1; ++j) {

      // The data frame j covers the range of phi such that
      // rj = {phi':phi0 + j*dphi <= phi' >= phi0 + (j+1)*dpi}
      // Therefore the range of phi for j is given as follows.
      double aj = starting_angle_ + j * oscillation_;
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
          double bv3 = ((v3 - 0.5) * step_size_e3_) / zeta + phi;
          double av3 = ((v3 + 0.5) * step_size_e3_) / zeta + phi;
          if (av3 > bv3) std::swap(av3, bv3);

          // We need to integrate over the intersection of sets rv3 and rj
          double av3j = std::max(av3, aj);
          double bv3j = std::min(bv3, bj);

          // If there is no intersection then set the fraction of the
          // counts contributed by data frame j to grid coordinate v3 to
          // zero, otherwise calculate it as the ratio of the integral
          // over the intersection or rv3 and rj to the integral over rj
          if (av3j >= bv3j) {
            fraction[i] = 0;
          } else {
            fraction[i] = (erf((bv3j - phi) * sigr2) -
                           erf((av3j - phi) * sigr2)) * integral_j_r;
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
   * A class to calculate the beam vector at each sub-divided pixel coordinate
   * that will be used during the transformation of the reflections from the
   * detector to the xds coordinate frame. This class is used during pre-
   * processing since no knowledge of the specific reflections are needed in
   * order to calculate the beam vectors. The beam vectors are then used along
   * with reflection specific stuff to calculate the xds coordinate for each
   * pixel.
   */
  class ReciprocalSpaceTransformDetectorLabCoords {

  public:

    /**
     * Initialise the class.
     */
    ReciprocalSpaceTransformDetectorLabCoords() {}

    /**
     * Calculate the beam vector at every pixel on the detector, sub-divided
     * into (n_div * n_div) equal areas. This is done to remove a certain
     * amount of processing from being done per reflection and ensuring it
     * is only done before the reflections are procesed.
     * @param detector The detector model
     * @param beam The beam model
     * @param n_div The number of sub-divisions to use
     * @returns An array of beam vectors
     */
    flex_vec3_double operator()(const Detector &detector,
                                const Beam &beam,
                                std::size_t n_div) const {

      // check the input
      DIALS_ASSERT(beam.get_wavelength() > 0.0);
      DIALS_ASSERT(n_div > 0);

      // Calculate the image size
      vec2<std::size_t> image_size = detector.get_image_size();
      std::size_t x_size = image_size[0] * n_div;
      std::size_t y_size = image_size[1] * n_div;
      double n_div_r = 1.0 / (double)n_div;
      double wavelength_r = 1.0 / beam.get_wavelength();

      // Create the necessary arrays
      flex_vec3_double detector_s1(flex_grid<>(y_size, x_size));

      // Calculate the beam vectors for each sub-division of the detector
      for (std::size_t j = 0, k = 0; j < y_size; ++j) {
        for (std::size_t i = 0; i < x_size; ++i, ++k) {
          detector_s1[k] = detector.get_pixel_lab_coord(
            vec2<double>((i + 0.5) * n_div_r,
                         (j + 0.5) * n_div_r)).normalize() * wavelength_r;
        }
      }

      // Return the s1 vector
      return detector_s1;
    }
  };

  /**
   * Class representing the XDS transform of the reflection profile on the
   * detector to the XDS reciprocal lattice coordinate frame.
   */
  class ReciprocalSpaceTransform {

  public:

    /**
     * Initialise the transform.
     * @param beam The beam model
     * @param detector The detector model
     * @param gonio The goniometer model
     * @param scan The scan model
     * @param mosaicity The crystal mosaicity
     * @param n_sigma The number of standard deviations to use
     * @param grid_half_size The size of the grid about the origin
     * @param n_div The number of pixel sub divisions to use (default 5)
     */
    ReciprocalSpaceTransform(const Beam &beam, const Detector &detector,
                             const Goniometer &gonio, const Scan &scan,
                             double mosaicity, double n_sigma,
                             std::size_t grid_half_size, std::size_t n_div = 5)
      : e3_fraction_(scan, mosaicity, n_sigma, grid_half_size),
        detector_s1_(ReciprocalSpaceTransformDetectorLabCoords()(
          detector, beam, n_div)),
        s0_(beam.get_s0()),
        m2_(gonio.get_rotation_axis()) {

      // Check some input
      DIALS_ASSERT(n_div > 0);
      DIALS_ASSERT(grid_half_size > 0);

      // Set the image volume size
      volume_ = vec3<std::size_t>(
          scan.get_array_range()[1],
          detector.get_image_size()[1],
          detector.get_image_size()[0]);

      // Set the grid size specifiers
      grid_half_size_ = grid_half_size;
      grid_size_ = 2 * grid_half_size + 1;

      // Calculate the grid step size
      double delta_d = beam.get_sigma_divergence() * n_sigma;
      double delta_m = mosaicity * n_sigma;
      step_size_ = vec3<double>(
        delta_d / grid_size_,
        delta_d / grid_size_,
        delta_m / grid_size_);

      // Check the step size is valid
      DIALS_ASSERT(step_size_[0] > 0 && step_size_[1] > 0 && step_size_[2] > 0);

      // Set the number of sub-divisions
      n_div_ = n_div;
    }

    flex_double operator()(const flex_int &pixels, const flex_int &mask,
                           int6 bbox, vec3<double> s1, double phi) const;

    /**
     * Transform a single reflection
     * @param reflection The reflection
     */
    void operator()(Reflection &reflection) const {

      // Get the shoebox profiles
      flex_int shoebox = reflection.get_shoebox();
      flex_int background = reflection.get_shoebox_background();
      flex_double transformed = reflection.get_transformed_shoebox();

      // Check they're the right size
      DIALS_ASSERT(transformed.accessor().all().all_eq(grid_size_));

      // Copy the background subtracted pixels to a new array
      flex_int shoebox_input(shoebox.accessor());
      for (std::size_t i = 0; i < shoebox.size(); ++i) {
        shoebox_input[i] = shoebox[i] - background[i];
      }

      // Transform the shoebox
      flex_double result = this->operator()(shoebox_input,
        reflection.get_shoebox_mask(), reflection.get_bounding_box(),
        reflection.get_beam_vector(), reflection.get_rotation_angle());

      // Copy the transformed shoebox back to the reflection
      for (std::size_t i = 0; i < transformed.size(); ++i) {
        transformed[i] = result[i];
      }
    }

    /**
     * Transform all the reflections in a list
     * @param reflections The reflection list
     */
    void operator()(ReflectionList &reflections) const {
      for (std::size_t i = 0; i < reflections.size(); ++i) {
        try {
          this->operator()(reflections[i]);
        } catch (dials::error) {
          continue;
        }
      }
    }

  private:

    ReciprocalSpaceTransformE3Fraction e3_fraction_;
    flex_vec3_double detector_s1_;
    vec3<double> s0_;
    vec3<double> m2_;
    vec3<std::size_t> volume_;
    std::size_t grid_half_size_;
    std::size_t grid_size_;
    vec3<double> step_size_;
    std::size_t n_div_;
  };

  /**
   * Transform the profile of the reflection at detector point xyz, with beam
   * vector s1 and rotation angle phi, to the XDS reciprocal lattice coordinate
   * frame.
   *
   * We treat the image pixels as bins of a histogram. In effect we want to
   * redistribute the image pixel counts to the XDS grid elements. We do this
   * by transforming the detector pixel coordinates around the reflection to the
   * XDS reciprocal lattice coordinate frame and determining the fraction of the
   * pixel value that is given to each element in the XDS grid.
   *
   * @param image The reflection image pixels
   * @param mask The reflection pixel mask
   * @param bbox The bounding box of the reflection
   * @param s1 The beam vector of the reflection
   * @param phi The rotation angle of the reflection
   * @returns The transformed grid
   * @throws std::rumtime_error if input is invalid
   */
  flex_double ReciprocalSpaceTransform::operator()(const flex_int &image,
      const flex_int &mask, int6 bbox, vec3<double> s1, double phi) const
  {
    // Constant for scaling values
    static const double r2d = 1.0 / scitbx::constants::pi_180;

    // Check the image sizes match
    DIALS_ASSERT(image.accessor().all().all_eq(mask.accessor().all()));

    // Check the bounding box size matches the image size
    DIALS_ASSERT(bbox[1] - bbox[0] == image.accessor().all()[2]);
    DIALS_ASSERT(bbox[3] - bbox[2] == image.accessor().all()[1]);
    DIALS_ASSERT(bbox[5] - bbox[4] == image.accessor().all()[0]);

    // Check the bounding box is in the image range
    DIALS_ASSERT(bbox[0] >= 0 && bbox[1] <= volume_[0]);
    DIALS_ASSERT(bbox[3] >= 0 && bbox[3] <= volume_[1]);
    DIALS_ASSERT(bbox[4] >= 0 && bbox[5] <= volume_[2]);

    // Get the grid data array
    flex_double grid(flex_grid<>(grid_size_, grid_size_, grid_size_));

    // Calculate the x, y, z ranges to iterate over
    std::size_t x0 = bbox[0] * n_div_, x1 = bbox[1] * n_div_;
    std::size_t y0 = bbox[2] * n_div_, y1 = bbox[3] * n_div_;
    std::size_t z0 = bbox[4], z1 = bbox[5];

    // Calculate 1 / n_div and 1 / (n_div*n_div) for convenience
    double n_div_r = 1.0 / n_div_;
    double div_fraction = n_div_r * n_div_r;

    // Calculate the reflection coordinate system e1 and e2 axes, and zeta, the
    // lorentz correction (used to calculate component on e3 axis
    vec3 <double> e1 = s1.cross(s0_).normalize();
    vec3 <double> e2 = s1.cross(e1).normalize();
    double zeta = m2_ * e1;
    double s1_length = s1.length();
    e1 = e1 * r2d / s1_length;
    e2 = e2 * r2d / s1_length;

    // Calculate the fraction of counts contributed by each data frame, j,
    // around the reflection to each grid point, v3 in the profile frame. Hold
    // these fractions in a 2d array.
    flex_double fraction = e3_fraction_(vec2<int>(z0, z1), phi, zeta);

    // Loop through all the pixels (and their sub-divisions). Calculate the
    // coordinate of each pixel in the XDS coordinate frame e1 and e2 axes.
    // Find the grid point in which the calculate point is contained and then
    // add the counts for that pixel to the grid. See Kabsch 2010
    for (std::size_t yy = y0; yy < y1; ++yy) {
      for (std::size_t xx = x0; xx < x1; ++xx) {
        vec3<double> ds = detector_s1_(yy, xx) - s1;
        double c1 = e1 * ds;
        double c2 = e2 * ds;
        int gi = grid_half_size_ + c1 / step_size_[2];
        int gj = grid_half_size_ + c2 / step_size_[1];
        if (gi < 0 || gi >= grid_size_ || gj < 0 || gj >= grid_size_) {
          continue;
        }
        std::size_t x = xx * n_div_r;
        std::size_t y = yy * n_div_r;
        for (std::size_t z = z0; z <= z1; ++z) {
          if (mask(z, y, z) != 0) {
            int value = image(z, y, x) * div_fraction;
            for (std::size_t gk = 0; gk < grid_size_; ++gk) {
              grid(gk, gj, gi) += value * fraction(z, gk);
            }
          }
        }
      }
    }

    // Return the grid
    return grid;
  }

}} // namespace = dials::algorithms

#endif // DIALS_ALGORITHMS_INTEGRATION_RECIPROCAL_SPACE_TRANSFORM_H
