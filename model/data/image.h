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

namespace dials { namespace model {

  /**
   * A class to pass multi panel images from python to c++
   */
  class Image {
  public:

    typedef af::versa< int, af::c_grid<2> > image_type;
    typedef af::ref< int, af::c_grid<2> > image_ref_type;
    typedef af::const_ref< int, af::c_grid<2> > image_const_ref_type;

    Image(image_type data)
      : data_(1) {
      data_[0] = data;
    }

    Image(const af::const_ref<image_type> &data)
      : data_(data.begin(), data.end()) {}

    std::size_t npanels() const {
      return data_.size();
    }

    image_const_ref_type operator[](std::size_t panel) const {
      return data_[panel].const_ref();
    }

  private:

    af::shared<image_type> data_;
  };

}}

#endif // DIALS_MODEL_DATA_IMAGE_H
