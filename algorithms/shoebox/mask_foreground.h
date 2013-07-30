/*
 * mask_foreground.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_SHOEBOX_MASK_FOREGROUND_H
#define DIALS_ALGORITHMS_SHOEBOX_MASK_FOREGROUND_H

namespace dials { namespace algorithms { namespace shoebox {

  /**
   * A class to mask foreground/background pixels
   */
  class MaskForeground {
  public:

    /**
     * Initialise the stuff needed to create the mask.
     * @param beam The beam model
     * @param detector The detector model
     * @param gonio The goniometer model
     * @param scan The scan model
     * @param delta_b nsigma * sigma_divergence
     * @param delta_m nsigma * mosaicity
     */
    MaskForeground(const Beam &beam, const Detector &detector,
                   const Goniometer &gonio, const Scan &scan,
                   double delta_b, double delta_m)
      : s1_map_(beam_vector_map(detector, beam, false)),
        m2_(gonio.get_rotation_axis()),
        s0_(beam.get_s0()),
        phi0_(scan.get_oscillation()[0]),
        dphi_(scan.get_oscillation()[1]),
        index0_(scan.get_array_range()[0]),
        delta_b_r_(1.0 / delta_b),
        delta_m_r_(1.0 / delta_m) {}

    /**
     * Set all the foreground/background pixels in the reflection mask.
     * @param reflection The reflection to mask
     */
    void operator()(Reflection &reflection) const {

      // Ensure the reflection is valid
      if (reflection.is_valid()) {

        // Get some bits from the reflection
        flex_int mask = reflection.get_shoebox_mask();
        vec3<double> s1 = reflection.get_beam_vector();
        double phi = reflection.get_rotation_angle();

        // Create the coordinate system and generators
        CoordinateSystem cs(m2_, s0_, s1, phi);
        CoordinateGenerator coordxy(cs, x0, y0, s1_map_);
        FromRotationAngleFast coordz(cs);

        // Loop through all the pixels in the shoebox, transform the point
        // to the reciprocal space coordinate system and check that it is
        // within the ellipse defined by:
        // (c1 / delta_b)^2 + (c2 / delta_b)^2 + (c3 / delta_m)^2 <= 1
        // Mark those points within as Foreground and those without as
        // Background.
        for (std::size_t j = 0; j < ysize; ++j) {
          for (std::size_t i = 0; i < xsize; ++i) {
            vec2<double> gxy = coordxy(j, i)
            double gxa2 = (gxy[0] * delta_b_r_)*(gxy[0] * delta_b_r_);
            double gyb2 = (gxy[1] * delta_b_r_)*(gxy[1] * delta_b_r_);
            for (std::size_t k = 0; k < zsize; ++k) {
              phi_dash = phi0_ + (k + z0 - index0_) * dphi_;
              double gz = coordz(phi_dash);
              double gzc2 = (gz * delta_m_r_)*(gz * delta_m_r_);
              if (gxa2 + gyb2 + gzc2 <= 1.0) {
                mask(k, j, i) |= Foreground;
              } else {
                mask(k, j, i) |= Background;
              }
            }
          }
        }
      }
    }

    /**
     * Mask all the foreground/background pixels for all the reflections
     * @param reflections The reflection list
     */
    void operator()(ReflectionList &reflections) const {
      for (std::size_t i = 0; i < reflections.size(); ++i) {
        this->operator()(reflections[i]);
      }
    }

  private:
    flex_vec3_double s1_map_;
    vec3<double> m2_;
    vec3<double> s0_;
    double phi0_;
    double dphi_;
    int index0_;
    double delta_b_r_;
    double delta_m_r_;
  };

}}} // namespace dials::algorithms::shoebox

#endif /* DIALS_ALGORITHMS_SHOEBOX_MASK_FOREGROUND_H */
