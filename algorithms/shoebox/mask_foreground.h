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

#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include <dxtbx/model/beam.h>
#include <dxtbx/model/goniometer.h>
#include <dxtbx/model/detector.h>
#include <dxtbx/model/scan.h>
#include <dials/model/data/shoebox.h>
#include <dials/algorithms/reflection_basis/coordinate_system.h>
#include <dials/algorithms/reflection_basis/beam_vector_map.h>
#include <dials/algorithms/shoebox/mask_code.h>
#include <dials/error.h>

namespace dials { namespace algorithms { namespace shoebox {

  using scitbx::vec2;
  using scitbx::vec3;
  using scitbx::af::int6;
  using dxtbx::model::Beam;
  using dxtbx::model::Detector;
  using dxtbx::model::Goniometer;
  using dxtbx::model::Scan;
  using dials::model::Shoebox;
  using dials::algorithms::reflection_basis::CoordinateSystem;
  using dials::algorithms::reflection_basis::transform::beam_vector_map;

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
                   const af::const_ref<double> &delta_b,
                   const af::const_ref<double> &delta_m)
      : m2_(gonio.get_rotation_axis()),
        s0_(beam.get_s0()),
        phi0_(scan.get_oscillation()[0]),
        dphi_(scan.get_oscillation()[1]),
        index0_(scan.get_array_range()[0]),
        index1_(scan.get_array_range()[1]){
      DIALS_ASSERT(delta_b.size() == delta_m.size());
      DIALS_ASSERT(delta_b.size() == scan.get_num_images());
      DIALS_ASSERT(delta_b.all_gt(0.0));
      DIALS_ASSERT(delta_m.all_gt(0.0));
      for (std::size_t i = 0; i < detector.size(); ++i) {
        s1_map_.push_back(beam_vector_map(detector[i], beam, true));
      }
      delta_b_r_.resize(delta_b.size());
      delta_m_r_.resize(delta_m.size());
      for (std::size_t i = 0; i < delta_b.size(); ++i) {
        delta_b_r_[i] = 1.0 / delta_b[i];
        delta_m_r_[i] = 1.0 / delta_m[i];
      }
    }

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
      : m2_(gonio.get_rotation_axis()),
        s0_(beam.get_s0()),
        phi0_(scan.get_oscillation()[0]),
        dphi_(scan.get_oscillation()[1]),
        index0_(scan.get_array_range()[0]),
        index1_(scan.get_array_range()[1]){
      DIALS_ASSERT(delta_b > 0.0);
      DIALS_ASSERT(delta_m > 0.0);
      for (std::size_t i = 0; i < detector.size(); ++i) {
        s1_map_.push_back(beam_vector_map(detector[i], beam, true));
      }
      delta_b_r_.resize(scan.get_num_images());
      delta_m_r_.resize(scan.get_num_images());
      for (std::size_t i = 0; i < scan.get_num_images(); ++i) {
        delta_b_r_[i] = 1.0 / delta_b;
        delta_m_r_[i] = 1.0 / delta_m;
      }
    }

    /**
     * Set all the foreground/background pixels in the shoebox mask.
     * @param shoebox The shoebox to mask
     * @param s1 The beam vector
     * @param frame The frame number
     * @param panel The panel number
     */
    void single(Shoebox<> &shoebox, vec3<double> s1, double frame, std::size_t panel) const {

      // Ensure the shoebox is valid
      //if (shoebox.is_valid()) {

        // Get some bits from the shoebox
        af::ref< int, af::c_grid<3> > mask = shoebox.mask.ref();
        int6 bbox = shoebox.bbox;
        double phi = phi0_ + (frame - index0_) * dphi_;
        int x0 = bbox[0], x1 = bbox[1];
        int y0 = bbox[2], y1 = bbox[3];
        int z0 = bbox[4], z1 = bbox[5];
        int xsize = x1 - x0;
        int ysize = y1 - y0;
        int zsize = z1 - z0;

        int z = (int)floor(frame);
        DIALS_ASSERT(z >= index0_);
        DIALS_ASSERT(z < index1_);
        DIALS_ASSERT(z >= z0 && z < z1);
        double delta_b_r2 = delta_b_r_[z-index0_] * delta_b_r_[z-index0_];
        double delta_m_r2 = delta_m_r_[z-index0_] * delta_m_r_[z-index0_];

        // Get the beam vector map
        DIALS_ASSERT(panel < s1_map_.size());
        af::const_ref< vec3<double>, af::c_grid<2> > s1_map = s1_map_[panel].const_ref();

        // Check the size of the mask
        DIALS_ASSERT(mask.accessor()[0] == zsize);
        DIALS_ASSERT(mask.accessor()[1] == ysize);
        DIALS_ASSERT(mask.accessor()[2] == xsize);

        // Create the coordinate system and generators
        CoordinateSystem cs(m2_, s0_, s1, phi);

        // Get the size of the image
        std::size_t width = s1_map.accessor()[1];
        std::size_t height = s1_map.accessor()[0];

        // Loop through all the pixels in the shoebox, transform the point
        // to the reciprocal space coordinate system and check that it is
        // within the ellipse defined by:
        // (c1 / delta_b)^2 + (c2 / delta_b)^2 <= 1
        // Mark those points within as Foreground and those without as
        // Background.
        for (int j = 0; j < ysize; ++j) {
          for (int i = 0; i < xsize; ++i) {
            if (x0 + i >= 0 && y0 + j >= 0 &&
                x0 + i < width && y0 + j < height) {
              vec2<double> gxy1 = cs.from_beam_vector(s1_map(y0 + j, x0 + i));
              vec2<double> gxy2 = cs.from_beam_vector(s1_map(y0 + j + 1, x0 + i));
              vec2<double> gxy3 = cs.from_beam_vector(s1_map(y0 + j, x0 + i + 1));
              vec2<double> gxy4 = cs.from_beam_vector(s1_map(y0 + j + 1, x0 + i + 1));
              double dxy1 = (gxy1[0]*gxy1[0] + gxy1[1]*gxy1[1]) * delta_b_r2;
              double dxy2 = (gxy2[0]*gxy2[0] + gxy2[1]*gxy2[1]) * delta_b_r2;
              double dxy3 = (gxy3[0]*gxy3[0] + gxy3[1]*gxy3[1]) * delta_b_r2;
              double dxy4 = (gxy4[0]*gxy4[0] + gxy4[1]*gxy4[1]) * delta_b_r2;
              double dxy = std::min(std::min(dxy1, dxy2), std::min(dxy3, dxy4));
              for (std::size_t k = 0; k < zsize; ++k) {
                if (z0 + k >= index0_ && z0 + k < index1_) {
                  double gz1 = cs.from_rotation_angle_fast(phi0_ + (z0 + k - index0_) * dphi_);
                  double gz2 = cs.from_rotation_angle_fast(phi0_ + (z0 + k + 1 - index0_) * dphi_);
                  double gz = std::abs(gz1) < std::abs(gz2) ? gz1 : gz2;
                  double gzc2 = gz*gz*delta_m_r2;
                  int mask_value = (dxy + gzc2 <= 1.0) ? Foreground : Background;
                  mask(k, j, i) |= mask_value;
                }
              }
            }
          }
        }
      //}
    }

    /**
     * Mask all the foreground/background pixels for all the shoeboxes
     * @param shoeboxes The shoebox list
     * @param s1 The list of beam vectors
     * @param frame The list of frame numbers
     * @param panel The list of panel numbers
     */
    void array(
        af::ref<Shoebox<> > shoeboxes,
        const af::const_ref< vec3<double> > &s1,
        const af::const_ref< double > &frame,
        const af::const_ref< std::size_t > &panel) const {
      DIALS_ASSERT(shoeboxes.size() == s1.size());
      DIALS_ASSERT(shoeboxes.size() == frame.size());
      DIALS_ASSERT(shoeboxes.size() == panel.size());
      for (std::size_t i = 0; i < shoeboxes.size(); ++i) {
        this->single(shoeboxes[i], s1[i], frame[i], panel[i]);
      }
    }

  private:
    std::vector< af::versa< vec3<double>, af::c_grid<2> > > s1_map_;
    vec3<double> m2_;
    vec3<double> s0_;
    double phi0_;
    double dphi_;
    int index0_;
    int index1_;
    af::shared<double> delta_b_r_;
    af::shared<double> delta_m_r_;
  };


  /**
   * Class to help compute bbox for multiple experiments.
   */
  class MaskMultiForeground {
  public:

    /**
     * Add a bbox calculator to the list.
     */
    void push_back(const MaskForeground &obj) {
      compute_.push_back(obj);
    }

    /**
     * Get the number of calculators.
     */
    std::size_t size() const {
      return compute_.size();
    }

    /**
     * Mask all the foreground/background pixels for all the shoeboxes
     * @param id The experiment id
     * @param shoeboxes The shoebox list
     * @param s1 The list of beam vectors
     * @param frame The list of frame numbers
     */
    void operator()(
        const af::const_ref<std::size_t> &id,
        af::ref<Shoebox<> > shoeboxes,
        const af::const_ref< vec3<double> > &s1,
        const af::const_ref< double > &frame,
        const af::const_ref< std::size_t > &panel) const {
      DIALS_ASSERT(shoeboxes.size() == id.size());
      DIALS_ASSERT(shoeboxes.size() == s1.size());
      DIALS_ASSERT(shoeboxes.size() == frame.size());
      DIALS_ASSERT(shoeboxes.size() == panel.size());
      for (std::size_t i = 0; i < shoeboxes.size(); ++i) {
        DIALS_ASSERT(id[i] < size());
        compute_[id[i]].single(shoeboxes[i], s1[i], frame[i], panel[i]);
      }
    }

  private:

    std::vector<MaskForeground> compute_;
  };

}}} // namespace dials::algorithms::shoebox

#endif /* DIALS_ALGORITHMS_SHOEBOX_MASK_FOREGROUND_H */
