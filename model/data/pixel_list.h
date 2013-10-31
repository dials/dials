/*
 * pixel_list.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_MODEL_DATA_PIXEL_LIST_H
#define DIALS_MODEL_DATA_PIXEL_LIST_H

#include <scitbx/vec3.h>
#include <scitbx/array_family/tiny_types.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/error.h>

namespace dials { namespace model {

  using scitbx::af::int2;
  using scitbx::vec3;

  class PixelList {
  public:
    PixelList() :
       size_(0, 0),
       first_frame_(0) {}

    PixelList(int2 size, int first_frame)
      : size_(size),
        first_frame_(first_frame),
        last_frame_(first_frame) {}

    PixelList(int2 size, int2 frame_range,
              af::shared<int> values,
              af::shared< vec3<int> > coords)
      : size_(size),
        first_frame_(frame_range[0]),
        last_frame_(frame_range[1]),
        coords_(coords),
        values_(values) {}

    /** @returns The size of the 2D image */
    int2 size() const {
      return size_;
    }

    std::size_t num_pixels() const {
      return values_.size();
    }

    int first_frame() const {
      return first_frame_;
    }

    int last_frame() const {
      return last_frame_;
    }

    /** @returns The frame range */
    int2 frame_range() const {
      DIALS_ASSERT(last_frame_ >= first_frame_);
      return int2(first_frame_, last_frame_);
    }

    /** @returns The number of frames */
    std::size_t num_frames() const {
      DIALS_ASSERT(last_frame_ >= first_frame_);
      return last_frame_ - first_frame_;
    }

    void add_image(const af::const_ref< int, af::c_grid<2> > &image,
                   const af::const_ref< bool, af::c_grid<2> > &mask) {
      DIALS_ASSERT(image.accessor().all_eq(mask.accessor()));
      DIALS_ASSERT(image.accessor().all_eq(size_));
      for (std::size_t j = 0; j < size_[0]; ++j) {
        for (std::size_t i = 0; i < size_[1]; ++i) {
          if (mask(j, i)) {
            coords_.push_back(vec3<int>(last_frame_, j, i));
            values_.push_back(image(j, i));
          }
        }
      }
      last_frame_++;
    }

    /**
     * @returns The list of valid point coordinates
     */
    af::shared< vec3<int> > coords() const {
      return coords_;
    }

    /**
     * @returns The list of valid point values
     */
    af::shared<int> values() const {
      return values_;
    }

  private:

    int2 size_;
    int first_frame_;
    int last_frame_;
    af::shared< vec3<int> > coords_;
    af::shared<int> values_;
  };

}} // namespace dials::model

#endif /* DIALS_MODEL_DATA_PIXEL_LIST_H */
