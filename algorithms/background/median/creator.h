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

#ifndef DIALS_ALGORITHMS_BACKGROUND_MEDIAN_CREATOR_H
#define DIALS_ALGORITHMS_BACKGROUND_MEDIAN_CREATOR_H

#include <dials/array_family/reflection_table.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/model/data/shoebox.h>
#include <dials/model/data/image_volume.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  using model::ImageVolume;
  using model::MultiPanelImageVolume;
  using model::Shoebox;

  namespace detail {

    template <typename T>
    T median(const af::const_ref<T> &x) {
      af::shared<T> temp(x.begin(), x.end());
      std::nth_element(temp.begin(), temp.begin() + temp.size() / 2, temp.end());
      return temp[temp.size() / 2];
    }

  }  // namespace detail

  template <typename T>
  void create_from_arrays(const af::const_ref<T, af::c_grid<3> > &data,
                          af::ref<T, af::c_grid<3> > background,
                          af::ref<int, af::c_grid<3> > mask) {
    DIALS_ASSERT(data.accessor().all_eq(background.accessor()));
    DIALS_ASSERT(data.accessor().all_eq(mask.accessor()));

    // Compute number of background pixels
    std::size_t num_background = 0;
    int mask_code = Valid | Background;
    for (std::size_t i = 0; i < mask.size(); ++i) {
      if ((mask[i] & mask_code) == mask_code) {
        num_background++;
      }
    }
    DIALS_ASSERT(num_background > 0);

    // Allocate some arrays
    af::shared<double> Y(num_background, 0);
    std::size_t j = 0;
    for (std::size_t i = 0; i < mask.size(); ++i) {
      if ((mask[i] & mask_code) == mask_code) {
        DIALS_ASSERT(j < Y.size());
        DIALS_ASSERT(data[i] >= 0);
        Y[j++] = data[i];
      }
    }
    DIALS_ASSERT(j == Y.size());

    // Compute the median value for the starting value
    double median = detail::median(Y.const_ref());

    // Fill in the background shoebox values
    for (std::size_t i = 0; i < data.size(); ++i) {
      background[i] = median;
      if ((mask[i] & mask_code) == mask_code) {
        mask[i] |= BackgroundUsed;
      }
    }
  }

  /**
   * Compute the background values
   * @param sbox The shoeboxes
   * @returns Success True/False
   */
  inline af::shared<bool> create_from_shoebox(af::ref<Shoebox<> > sbox) {
    af::shared<bool> success(sbox.size(), true);
    for (std::size_t i = 0; i < sbox.size(); ++i) {
      try {
        DIALS_ASSERT(sbox[i].is_consistent());
        create_from_arrays(
          sbox[i].data.const_ref(), sbox[i].background.ref(), sbox[i].mask.ref());
      } catch (scitbx::error) {
        success[i] = false;
      } catch (dials::error) {
        success[i] = false;
      }
    }
    return success;
  }

  /**
   * Compute the background values
   * @param reflections The reflection table
   * @param volume The image volume
   * @returns Success True/False
   */
  inline af::shared<bool> create_from_image_volume(af::reflection_table reflections,
                                                   MultiPanelImageVolume<> volume) {
    typedef MultiPanelImageVolume<>::float_type FloatType;

    DIALS_ASSERT(reflections.contains("bbox"));
    DIALS_ASSERT(reflections.contains("panel"));
    af::const_ref<int6> bbox = reflections["bbox"];
    af::const_ref<std::size_t> panel = reflections["panel"];
    af::shared<bool> success(bbox.size(), true);
    for (std::size_t i = 0; i < bbox.size(); ++i) {
      // Get the image volume
      ImageVolume<> v = volume.get(panel[i]);

      // Trim the bbox
      int6 b = v.trim_bbox(bbox[i]);

      // Extract from image volume
      af::versa<FloatType, af::c_grid<3> > data = v.extract_data(b);
      af::versa<FloatType, af::c_grid<3> > bgrd = v.extract_background(b);
      af::versa<int, af::c_grid<3> > mask = v.extract_mask(b, i);

      // Compute the background
      try {
        create_from_arrays(data.const_ref(), bgrd.ref(), mask.ref());

        // Need to set the background in volume
        v.set_background(b, bgrd.const_ref());
      } catch (scitbx::error) {
        success[i] = false;
      } catch (dials::error) {
        success[i] = false;
      }
    }
    return success;
  }

}}  // namespace dials::algorithms

#endif  // DIALS_ALGORITHMS_BACKGROUND_MEDIAN_CREATOR_H
