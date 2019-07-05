/*
 * image.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */

#ifndef DIALS_MODEL_DATA_IMAGE_H
#define DIALS_MODEL_DATA_IMAGE_H

#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/error.h>

namespace dials { namespace model {

  /**
   * A class to pass multi panel images from python to c++
   */
  template <typename T>
  class Image {
  public:
    typedef T value_type;

    typedef af::versa<bool, af::c_grid<2> > bool_type;
    typedef af::ref<bool, af::c_grid<2> > bool_ref_type;
    typedef af::const_ref<bool, af::c_grid<2> > bool_const_ref_type;

    typedef af::versa<T, af::c_grid<2> > data_type;
    typedef af::ref<T, af::c_grid<2> > data_ref_type;
    typedef af::const_ref<T, af::c_grid<2> > data_const_ref_type;

    Image(data_type data, bool_type mask) : data_(1), mask_(1) {
      data_[0] = data;
      mask_[0] = mask;
      DIALS_ASSERT(data.accessor().all_eq(mask.accessor()));
    }

    Image(const af::const_ref<data_type> &data, const af::const_ref<bool_type> &mask)
        : data_(data.begin(), data.end()), mask_(mask.begin(), mask.end()) {
      DIALS_ASSERT(data_.size() == mask_.size());
      for (std::size_t i = 0; i < data_.size(); ++i) {
        DIALS_ASSERT(data_[i].accessor().all_eq(mask_[i].accessor()));
      }
    }

    std::size_t npanels() const {
      return data_.size();
    }

    data_const_ref_type data(std::size_t panel) const {
      return data_[panel].const_ref();
    }

    bool_const_ref_type mask(std::size_t panel) const {
      return mask_[panel].const_ref();
    }

  private:
    af::shared<data_type> data_;
    af::shared<bool_type> mask_;
  };

}}  // namespace dials::model

#endif  // DIALS_MODEL_DATA_IMAGE_H
