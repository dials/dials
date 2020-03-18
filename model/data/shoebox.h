/*
 * shoebox.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_MODEL_DATA_SHOEBOX_H
#define DIALS_MODEL_DATA_SHOEBOX_H

#include <scitbx/array_family/tiny_types.h>
#include <scitbx/array_family/small.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/model/data/observation.h>
#include <dials/algorithms/image/centroid/centroid_image.h>
#include <dials/algorithms/image/centroid/centroid_masked_image.h>
#include <dials/algorithms/integration/sum/summation.h>
#include <dials/algorithms/integration/bayes/bayesian_integrator.h>
#include <dials/model/data/mask_code.h>
#include <dials/config.h>
#include <dials/error.h>

namespace dials { namespace model {

  using dials::algorithms::BayesianIntegrator;
  using dials::algorithms::CentroidImage3d;
  using dials::algorithms::CentroidMaskedImage3d;
  using dials::algorithms::Summation;
  using dials::model::Centroid;
  using dials::model::Foreground;
  using dials::model::Intensity;
  using dials::model::Valid;
  using scitbx::af::int2;
  using scitbx::af::int3;
  using scitbx::af::int6;
  using scitbx::af::small;

  /**
   * Helper function to take an image centroid and create something we
   * can return as the centroid
   */
  template <typename AlgorithmType>
  Centroid extract_centroid_object(const AlgorithmType &algorithm,
                                   const vec3<double> &ioffset) {
    // Try to get the variance and add 0.5 to the standard error
    Centroid result;
    result.px.position = algorithm.mean() + ioffset;
    try {
      result.px.variance = algorithm.unbiased_variance();
      result.px.std_err_sq = algorithm.mean_sq_error();
    } catch (dials::error) {
      result.px.variance = vec3<double>(0.0, 0.0, 0.0);
      result.px.std_err_sq = vec3<double>(1.0 / 12.0, 1.0 / 12.0, 1.0 / 12.0);
    }
    return result;
  }

  /**
   * A class to hold shoebox information
   */
  template <typename FloatType = ProfileFloatType>
  struct Shoebox {
    typedef FloatType float_type;

    std::size_t panel;                                ///< The detector panel
    int6 bbox;                                        ///< The bounding box
    bool flat;                                        ///< Is the shoebox flat
    af::versa<FloatType, af::c_grid<3> > data;        ///< The shoebox data
    af::versa<int, af::c_grid<3> > mask;              ///< The shoebox mask
    af::versa<FloatType, af::c_grid<3> > background;  ///< The shoebox background

    /**
     * Initialise the shoebox
     */
    Shoebox()
        : panel(0),
          bbox(0, 0, 0, 0, 0, 0),
          flat(false),
          data(af::c_grid<3>(0, 0, 0)),
          mask(af::c_grid<3>(0, 0, 0)),
          background(af::c_grid<3>(0, 0, 0)) {}

    /**
     * Initialise the shoebox
     * @param bbox_ The bounding box to initialise with
     */
    Shoebox(const int6 &bbox_)
        : panel(0),
          bbox(bbox_),
          flat(false),
          data(af::c_grid<3>(0, 0, 0)),
          mask(af::c_grid<3>(0, 0, 0)),
          background(af::c_grid<3>(0, 0, 0)) {}

    /**
     * Initialise the shoebox
     * @param panel_ The panel number
     * @param bbox_ The bounding box to initialise with
     */
    Shoebox(std::size_t panel_, const int6 &bbox_)
        : panel(panel_),
          bbox(bbox_),
          flat(false),
          data(af::c_grid<3>(0, 0, 0)),
          mask(af::c_grid<3>(0, 0, 0)),
          background(af::c_grid<3>(0, 0, 0)) {}

    /**
     * Initialise the shoebox
     * @param panel_ The panel number
     * @param bbox_ The bounding box to initialise with
     * @param flat_ The shoebix is flat
     */
    Shoebox(std::size_t panel_, const int6 &bbox_, bool flat_)
        : panel(panel_),
          bbox(bbox_),
          flat(flat_),
          data(af::c_grid<3>(0, 0, 0)),
          mask(af::c_grid<3>(0, 0, 0)),
          background(af::c_grid<3>(0, 0, 0)) {}

    /**
     * Allocate the mask and data from the bounding box
     */
    void allocate_with_value(int maskcode) {
      std::size_t zs = flat ? 1 : zsize();
      af::c_grid<3> accessor(zs, ysize(), xsize());
      data = af::versa<FloatType, af::c_grid<3> >(accessor, 0.0);
      mask = af::versa<int, af::c_grid<3> >(accessor, maskcode);
      background = af::versa<FloatType, af::c_grid<3> >(accessor, 0.0);
    }

    /**
     * Allocate the mask and data with mask code valid.
     */
    void allocate() {
      allocate_with_value(0);
    }

    /**
     * Deallocate the mask and data arrays
     */
    void deallocate() {
      af::c_grid<3> accessor(0, 0, 0);
      data = af::versa<FloatType, af::c_grid<3> >(accessor);
      mask = af::versa<int, af::c_grid<3> >(accessor);
      background = af::versa<FloatType, af::c_grid<3> >(accessor);
    }

    /** @returns The x offset */
    int xoffset() const {
      return bbox[0];
    }

    /** @returns The y offset */
    int yoffset() const {
      return bbox[2];
    }

    /** @returns The z offset */
    int zoffset() const {
      return bbox[4];
    }

    /** @returns The x size */
    std::size_t xsize() const {
      DIALS_ASSERT(bbox[1] >= bbox[0]);
      return (std::size_t)(bbox[1] - bbox[0]);
    }

    /** @returns The y size */
    std::size_t ysize() const {
      DIALS_ASSERT(bbox[3] >= bbox[2]);
      return (std::size_t)(bbox[3] - bbox[2]);
    }

    /** @returns The z size */
    std::size_t zsize() const {
      DIALS_ASSERT(bbox[5] >= bbox[4]);
      return (std::size_t)(bbox[5] - bbox[4]);
    }

    /** @returns The offset */
    int3 offset() const {
      return int3(zoffset(), yoffset(), xoffset());
    }

    /** @returns The size */
    int3 size() const {
      return int3(zsize(), ysize(), xsize());
    }

    /** @returns The flat size */
    int3 size_flat() const {
      return int3(1, ysize(), xsize());
    }

    /** @return True/False whether the array and bbox sizes are consistent */
    bool is_consistent() const {
      bool result = true;
      if (flat) {
        result = result && (data.accessor().all_eq(size_flat()));
        result = result && (mask.accessor().all_eq(size_flat()));
        result = result && (background.accessor().all_eq(size_flat()));
      } else {
        result = result && (data.accessor().all_eq(size()));
        result = result && (mask.accessor().all_eq(size()));
        result = result && (background.accessor().all_eq(size()));
      }
      return result;
    }

    /** @return True/False whether the shoeboxes are allocated */
    bool is_allocated() const {
      DIALS_ASSERT(data.accessor().all_eq(mask.accessor()));
      DIALS_ASSERT(data.accessor().all_eq(background.accessor()));
      if (data.size() == 0) {
        return false;
      } else {
        if (flat) {
          DIALS_ASSERT(data.accessor().all_eq(size_flat()));
        } else {
          DIALS_ASSERT(data.accessor().all_eq(size()));
        }
      }
      return true;
    }

    /**
     * Check if the bounding box has points outside the image range.
     * @param image_size The image size
     * @param scan_range The scan range
     * @returns True/False
     */
    bool is_bbox_within_image_volume(int2 image_size, int2 scan_range) const {
      return bbox[0] >= 0 && bbox[1] < image_size[1] && bbox[2] >= 0
             && bbox[3] < image_size[0] && bbox[4] >= scan_range[0]
             && bbox[5] < scan_range[1];
    }

    /**
     * Check if the bounding box has points that cover bad pixels
     * @param mask The mask array
     * @returns True/False
     */
    bool does_bbox_contain_bad_pixels(
      const af::const_ref<bool, af::c_grid<2> > &mask) const {
      std::size_t ysize = mask.accessor()[0];
      std::size_t xsize = mask.accessor()[1];
      int j0 = bbox[2] > 0 ? bbox[2] : 0;
      int j1 = bbox[3] < ysize ? bbox[3] : ysize;
      int i0 = bbox[0] > 0 ? bbox[0] : 0;
      int i1 = bbox[1] < xsize ? bbox[1] : xsize;
      for (int j = j0; j < j1; ++j) {
        for (int i = i0; i < i1; ++i) {
          if (mask(j, i) == false) {
            return true;
          }
        }
      }
      return false;
    }

    /**
     * Count the number of mask pixels with the given value
     * @param code The code
     * @returns The number of pixels with that code
     */
    int count_mask_values(int code) const {
      int count = 0;
      for (std::size_t i = 0; i < mask.size(); ++i) {
        if ((mask[i] & code) == code) {
          count++;
        }
      }
      return count;
    }

    /**
     * returns true if all the foreground pixels are valid
     * @returns a bool
     */
    bool all_foreground_valid() const {
      for (std::size_t i = 0; i < mask.size(); ++i) {
        if (((mask[i] & Valid) != Valid) && ((mask[i] & Foreground) == Foreground)) {
          return false;
        }
      }
      return true;
    }

    /**
     * Perform a centroid of all pixels
     * @returns The centroid
     */
    Centroid centroid_all() const {
      typedef CentroidImage3d<FloatType> Centroider;
      Centroider centroid(data.const_ref());
      vec3<double> offset(bbox[0], bbox[2], bbox[4]);
      return extract_centroid_object(centroid, offset);
    }

    /**
     * Perform a centroid of masked pixels
     * @param code The mask code
     * @returns The centroid
     */
    Centroid centroid_masked(int code) const {
      typedef CentroidMaskedImage3d<FloatType> Centroider;

      // Calculate the foreground mask
      af::versa<bool, af::c_grid<3> > fg_mask_arr(mask.accessor());
      af::ref<bool, af::c_grid<3> > foreground_mask = fg_mask_arr.ref();
      for (std::size_t i = 0; i < mask.size(); ++i) {
        foreground_mask[i] = (mask[i] & code) == code && ((mask[i] & Overlapped) == 0);
      }

      // Calculate the centroid
      int zoff = flat ? (bbox[5] + bbox[4]) / 2 : bbox[4];
      vec3<double> offset(bbox[0], bbox[2], zoff);
      Centroid result;
      try {
        Centroider centroid(data.const_ref(), foreground_mask);
        result = extract_centroid_object(centroid, offset);
        if (bbox[5] == bbox[4] + 1) {
          result.px.position[2] = bbox[4] + 0.5;
        }
      } catch (dials::error) {
        double xmid = (bbox[1] + bbox[0]) / 2.0;
        double ymid = (bbox[3] + bbox[2]) / 2.0;
        double zmid = (bbox[5] + bbox[4]) / 2.0;
        result.px.position = vec3<double>(xmid, ymid, zmid);
        result.px.variance = vec3<double>(0, 0, 0);
        result.px.std_err_sq = vec3<double>(0, 0, 0);
      }
      return result;
    }

    /**
     * Perform a centroid of valid pixels
     * @returns The centroid
     */
    Centroid centroid_valid() const {
      return centroid_masked(Valid);
    }

    /**
     * Perform a centroid of foreground pixels
     * @returns The centroid
     */
    Centroid centroid_foreground() const {
      return centroid_masked(Valid | Foreground);
    }

    /**
     * Perform a centroid of strong pixels
     * @returns The centroid
     */
    Centroid centroid_strong() const {
      return centroid_masked(Valid | Strong);
    }

    /**
     * Perform a centroid minus the background
     * @return The centroid
     */
    Centroid centroid_all_minus_background() const {
      typedef CentroidImage3d<FloatType> Centroider;

      // Calculate the foreground mask and data
      DIALS_ASSERT(data.size() == background.size());
      af::versa<FloatType, af::c_grid<3> > fg_data_arr(data.accessor());
      af::ref<FloatType, af::c_grid<3> > foreground_data = fg_data_arr.ref();
      for (std::size_t i = 0; i < mask.size(); ++i) {
        foreground_data[i] = data[i] - background[i];
      }

      // Calculate the centroid
      Centroider centroid(foreground_data);
      int zoff = flat ? (bbox[5] + bbox[4]) / 2 : bbox[4];
      vec3<double> offset(bbox[0], bbox[2], zoff);
      return extract_centroid_object(centroid, offset);
    }

    /**
     * Perform a centroid minus the background
     * @return The centroid
     */
    Centroid centroid_masked_minus_background(int code) const {
      typedef CentroidMaskedImage3d<FloatType> Centroider;

      // Calculate the foreground mask and data
      DIALS_ASSERT(data.size() == mask.size());
      DIALS_ASSERT(data.size() == background.size());
      af::versa<bool, af::c_grid<3> > fg_mask_arr(mask.accessor());
      af::versa<FloatType, af::c_grid<3> > fg_data_arr(data.accessor());
      af::ref<FloatType, af::c_grid<3> > foreground_data = fg_data_arr.ref();
      af::ref<bool, af::c_grid<3> > foreground_mask = fg_mask_arr.ref();
      for (std::size_t i = 0; i < mask.size(); ++i) {
        foreground_data[i] = data[i] - background[i];
        foreground_mask[i] = ((mask[i] & code) == code) && ((mask[i] & Overlapped) == 0)
                             && (foreground_data[i] > 0);
      }

      // Calculate the centroid
      int zoff = flat ? (bbox[5] + bbox[4]) / 2 : bbox[4];
      vec3<double> offset(bbox[0], bbox[2], zoff);
      Centroid result;
      try {
        Centroider centroid(foreground_data, foreground_mask);
        result = extract_centroid_object(centroid, offset);
      } catch (dials::error) {
        double xmid = (bbox[1] + bbox[0]) / 2.0;
        double ymid = (bbox[3] + bbox[2]) / 2.0;
        double zmid = (bbox[5] + bbox[4]) / 2.0;
        result.px.position = vec3<double>(xmid, ymid, zmid);
        result.px.variance = vec3<double>(0, 0, 0);
        result.px.std_err_sq = vec3<double>(0, 0, 0);
      }
      return result;
    }

    /**
     * Perform a centroid minus the background
     * @return The centroid
     */
    Centroid centroid_valid_minus_background() const {
      return centroid_masked_minus_background(Valid);
    }

    /**
     * Perform a centroid minus the background
     * @return The centroid
     */
    Centroid centroid_foreground_minus_background() const {
      return centroid_masked_minus_background(Valid | Foreground);
    }

    /**
     * Perform a centroid minus the background
     * @return The centroid
     */
    Centroid centroid_strong_minus_background() const {
      return centroid_masked_minus_background(Valid | Strong);
    }

    /**
     * Get the summed intensity of all pixels
     * @returns The intensity
     */
    Intensity bayesian_intensity() const {
      // Do the intengration
      BayesianIntegrator<FloatType> summation(
        data.const_ref(), background.const_ref(), mask.const_ref());

      // Return the intensity struct
      Intensity result;
      result.observed.value = summation.intensity();
      result.observed.variance = summation.variance();
      result.background.value = summation.background();
      result.background.variance = summation.background_variance();
      result.observed.success = summation.success();
      return result;
    }

    /**
     * Get the summed intensity of all pixels
     * @returns The intensity
     */
    Intensity summed_intensity() const {
      // Do the intengration
      Summation<FloatType> summation(
        data.const_ref(), background.const_ref(), mask.const_ref());

      // Return the intensity struct
      Intensity result;
      result.observed.value = summation.intensity();
      result.observed.variance = summation.variance();
      result.background.value = summation.background();
      result.background.variance = summation.background_variance();
      result.observed.success = summation.success();
      return result;
    }

    /**
     * Test to see if shoeboxs contain the same data
     * @param rhs The other shoebox
     * @returns True/False. They are the same
     */
    bool operator==(const Shoebox &rhs) const {
      return ((bbox.all_eq(rhs.bbox)) && (data.all_eq(rhs.data))
              && (mask.all_eq(rhs.mask)) && (background.all_eq(rhs.background)));
    }

    /**
     * Test to see if shoeboxs contain the same data
     * @param rhs The other shoebox
     * @returns True/False. They are not the same
     */
    bool operator!=(const Shoebox &rhs) const {
      return !(*this == rhs);
    }

    /**
     * Flatten the shoebox
     */
    void flatten() {
      DIALS_ASSERT(is_consistent());
      if (flat == false) {
        for (std::size_t k = 1; k < zsize(); ++k) {
          for (std::size_t j = 0; j < ysize(); ++j) {
            for (std::size_t i = 0; i < xsize(); ++i) {
              data(0, j, i) += data(k, j, i);
              mask(0, j, i) |= mask(k, j, i);
              bool valid = (mask(0, j, i) & Valid) && (mask(0, j, i) & Valid);
              if (!valid) {
                mask(0, j, i) &= ~Valid;
              }
            }
          }
        }
        af::c_grid<3> accessor(1, ysize(), xsize());
        data.resize(accessor);
        mask.resize(accessor);
        background.resize(accessor);
      }
      flat = true;
      DIALS_ASSERT(is_consistent());
    }
  };

}};  // namespace dials::model

#endif /* DIALS_MODEL_DATA_SHOEBOX_H */
