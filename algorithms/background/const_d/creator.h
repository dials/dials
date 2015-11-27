/*
 * creator.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */

#ifndef DIALS_ALGORITHMS_BACKGROUND_GLM_CREATOR_H
#define DIALS_ALGORITHMS_BACKGROUND_GLM_CREATOR_H

#include <dxtbx/model/detector.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/model/data/shoebox.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  using dxtbx::model::Beam;
  using dxtbx::model::Detector;
  using dxtbx::model::Panel;
  using model::Shoebox;

  template <typename T>
  T min4(T a, T b, T c, T d) {
    return std::min(std::min(a,b), std::min(c,d));
  }

  template <typename T>
  T max4(T a, T b, T c, T d) {
    return std::max(std::max(a,b), std::max(c,d));
  }

  /**
   * A class to create the background model
   */
  class Creator {
  public:

    /**
     * Initialise the creator
     */
    Creator(const Beam &beam, const Detector &detector)
      : beam_(beam),
        detector_(detector) {}

    /**
     * Compute the background values
     * @param sbox The shoeboxes
     * @returns Success True/False
     */
    af::shared<bool> operator()(af::ref< Shoebox<> > sbox) const {
      af::shared<bool> success(sbox.size(), true);
      for (std::size_t i = 0; i < sbox.size(); ++i) {
        try {
          compute(sbox[i]);
        } catch(scitbx::error) {
          success[i] = false;
        } catch(dials::error) {
          success[i] = false;
        }
      }
      return success;
    }

  private:

    /**
     * Compute the background values for a single shoebox
     * @param sbox The shoebox
     */
    void compute(Shoebox<> &sbox) const {
      DIALS_ASSERT(sbox.is_consistent());

      // Get the panel
      Panel panel = detector_[sbox.panel];
      vec3<double> s0 = beam_.get_s0();

      // Compute the resolution at each pixel
      std::vector<double> d_min(sbox.data.size());
      std::vector<double> d_max(sbox.data.size());
      std::vector<double> d_mid(sbox.data.size());
      for (std::size_t k = 0; k < sbox.zsize(); ++k) {
        for (std::size_t j = 0; j < sbox.ysize(); ++j) {
          for (std::size_t i = 0; i < sbox.xsize(); ++i) {
            double x = i + sbox.bbox[0];
            double y = j + sbox.bbox[2];
            vec2<double> xy1(x,y);
            vec2<double> xy2(x+1,y);
            vec2<double> xy3(x,y+1);
            vec2<double> xy4(x+1,y+1);
            vec2<double> xy5(x+0.5,y+0.5);
            double d1 = panel.get_resolution_at_pixel(s0, xy1);
            double d2 = panel.get_resolution_at_pixel(s0, xy2);
            double d3 = panel.get_resolution_at_pixel(s0, xy3);
            double d4 = panel.get_resolution_at_pixel(s0, xy4);
            double d5 = panel.get_resolution_at_pixel(s0, xy5);
            d_min.push_back(min4(d1,d2,d3,d4));
            d_max.push_back(max4(d1,d2,d3,d4));
            d_mid.push_back(d5);
          }
        }
      }

      // Get the indices of background and foreground pixels
      int background_code = Valid | Background;
      int foreground_code = Valid | Foreground;
      std::vector<std::size_t> background;
      std::vector<std::size_t> foreground;
      af::const_ref< int, af::c_grid<3> > mask = sbox.mask.const_ref();
      for (std::size_t i = 0; i < mask.size(); ++i) {
        if ((mask[i] & background_code) == background_code) {
          background.push_back(i);
        }
        if ((mask[i] & foreground_code) == foreground_code) {
          foreground.push_back(i);
        }
      }

      // Compute radially averaged background
      af::const_ref< float, af::c_grid<3> > data = sbox.data.const_ref();
      af::ref< float, af::c_grid<3> > background_values = sbox.data.ref();
      for (std::size_t i = 0; i < foreground.size(); ++i) {
        std::size_t f = foreground[i];
        double dm = d_mid[f];
        double sum = 0.0;
        std::size_t count = 0;
        for (std::size_t j = 0; j < background.size(); ++j) {
          std::size_t b = background[i];
          double d0 = d_min[b];
          double d1 = d_max[b];
          if (d0 <= dm && dm <= d1) {
            sum += data[b];
            count += 1;
          }
        }
        if (count > 0) {
          background_values[f] = sum / count;
        } else {
          background_values[f] = 0;
        }
      }
    }

    Beam beam_;
    Detector detector_;
  };

}} // namespace dials::algorithms

#endif // DIALS_ALGORITHMS_BACKGROUND_GLM_CREATOR_H
