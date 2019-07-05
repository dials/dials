/*
 * centroid.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */

#ifndef DIALS_ALGORITHMS_CENTROID_CENTROID_H
#define DIALS_ALGORITHMS_CENTROID_CENTROID_H

#include <dials/algorithms/image/centroid/centroid_masked_image.h>
#include <dials/model/data/shoebox.h>

namespace dials { namespace algorithms {

  using model::Foreground;
  using model::Valid;

  class Centroid {
  public:
    Centroid(const af::const_ref<double, af::c_grid<2> > &data,
             const af::const_ref<bool, af::c_grid<2> > &mask) {
      DIALS_ASSERT(data.accessor().all_eq(mask.accessor()));
      init(data, mask);
    }

    Centroid(const af::const_ref<double, af::c_grid<3> > &data,
             const af::const_ref<bool, af::c_grid<3> > &mask) {
      DIALS_ASSERT(data.accessor().all_eq(mask.accessor()));
      init(data, mask);
    }

    Centroid(const af::const_ref<double, af::c_grid<2> > &data,
             const af::const_ref<double, af::c_grid<2> > &background,
             const af::const_ref<int, af::c_grid<2> > &mask) {
      DIALS_ASSERT(data.accessor().all_eq(background.accessor()));
      DIALS_ASSERT(data.accessor().all_eq(mask.accessor()));
      init(data, background, mask);
    }

    Centroid(const af::const_ref<double, af::c_grid<3> > &data,
             const af::const_ref<double, af::c_grid<3> > &background,
             const af::const_ref<int, af::c_grid<3> > &mask) {
      DIALS_ASSERT(data.accessor().all_eq(background.accessor()));
      DIALS_ASSERT(data.accessor().all_eq(mask.accessor()));
      init(data, background, mask);
    }

    vec3<double> value() const {
      return value_;
    }

    vec3<double> variance() const {
      return variance_;
    }

  private:
    void init(const af::const_ref<double, af::c_grid<2> > &data,
              const af::const_ref<bool, af::c_grid<2> > &mask) {
      CentroidMaskedImage2d<> centroider(data, mask);
      vec2<double> m = centroider.mean();
      vec2<double> v = centroider.unbiased_standard_error_sq();
      value_[0] = m[0];
      value_[1] = m[1];
      value_[2] = 0;
      variance_[0] = v[0];
      variance_[1] = v[1];
      variance_[2] = 0;
    }

    void init(const af::const_ref<double, af::c_grid<3> > &data,
              const af::const_ref<bool, af::c_grid<3> > &mask) {
      CentroidMaskedImage3d<> centroider(data, mask);
      value_ = centroider.mean();
      variance_ = centroider.unbiased_standard_error_sq();
    }

    void init(const af::const_ref<double, af::c_grid<2> > &data,
              const af::const_ref<double, af::c_grid<2> > &background,
              const af::const_ref<int, af::c_grid<2> > &mask) {
      af::versa<double, af::c_grid<2> > temp_data(data.accessor());
      af::versa<bool, af::c_grid<2> > temp_mask(data.accessor());
      int fg_code = Valid | Foreground;
      for (std::size_t i = 0; i < temp_data.size(); ++i) {
        temp_data[i] = data[i] - background[i];
        temp_mask[i] = (mask[i] & fg_code) == fg_code ? true : false;
      }
      init(temp_data.const_ref(), temp_mask.const_ref());
    }

    void init(const af::const_ref<double, af::c_grid<3> > &data,
              const af::const_ref<double, af::c_grid<3> > &background,
              const af::const_ref<int, af::c_grid<3> > &mask) {
      af::versa<double, af::c_grid<3> > temp_data(data.accessor());
      af::versa<bool, af::c_grid<3> > temp_mask(data.accessor());
      int fg_code = Valid | Foreground;
      for (std::size_t i = 0; i < temp_data.size(); ++i) {
        temp_data[i] = data[i] - background[i];
        temp_mask[i] = (mask[i] & fg_code) == fg_code ? true : false;
      }
      init(temp_data.const_ref(), temp_mask.const_ref());
    }

    vec3<double> value_;
    vec3<double> variance_;
  };

}}  // namespace dials::algorithms

#endif  // DIALS_ALGORITHMS_CENTROID_CENTROID_H
