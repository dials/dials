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
#include <dials/algorithms/image/connected_components/connected_components.h>
#include <dials/error.h>

namespace dials { namespace model {

  using scitbx::af::int2;
  using scitbx::vec3;
  using dials::algorithms::labelpixels3d;
  using dials::algorithms::labelpixels2d;

  /**
   * A class to hold a load of pixels and coordinates
   */
  class PixelList {
  public:

    PixelList() :
       size_(0, 0),
       first_frame_(0) {}

    /**
     * Initialise
     * @param size The size of the 2D image
     * @param first_frame The index of the first frame
     */
    PixelList(int2 size, int first_frame)
      : size_(size),
        first_frame_(first_frame),
        last_frame_(first_frame) {}

    /**
     * Initialise
     * @param size The size of the 2D image
     * @param frame_range The range of the frame indices
     * @param values The image values
     * @param coords The image coords
     */
    PixelList(int2 size, int2 frame_range,
              af::shared<int> values,
              af::shared< vec3<int> > coords)
      : size_(size),
        first_frame_(frame_range[0]),
        last_frame_(frame_range[1]),
        coords_(coords),
        values_(values) {
      DIALS_ASSERT(last_frame_ >= first_frame_);
      DIALS_ASSERT(coords_.size() == values_.size());
      for (std::size_t i = 0; i < coords_.size(); ++i) {
        DIALS_ASSERT(coords[i][0] >= first_frame_ && coords[i][0] < last_frame_);
        DIALS_ASSERT(coords[i][1] >= 0 && coords[i][1] < size_[0]);
        DIALS_ASSERT(coords[i][2] >= 0 && coords[i][2] < size_[1]);
      }
    }

    /** @returns The size of the 2D image */
    int2 size() const {
      return size_;
    }

    /** @returns The number of pixels */
    std::size_t num_pixels() const {
      return values_.size();
    }

    /** @returns The first frame number */
    int first_frame() const {
      return first_frame_;
    }

    /** @returns The last frame number */
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

    /**
     * Add an image
     * @param image The image pixels
     * @param mask The mask values
     */
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

    af::shared<int> labels_3d() const {
      return labelpixels3d(coords_.const_ref());
    }

    af::shared<int> labels_2d() const {
      return labelpixels2d(coords_.const_ref());
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
