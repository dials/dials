/*
 * image_volume.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_MODEL_DATA_IMAGE_VOLUME_H
#define DIALS_MODEL_DATA_IMAGE_VOLUME_H

#include <boost/unordered_map.hpp>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/model/data/image.h>
#include <dials/model/data/mask_code.h>
#include <dials/error.h>

namespace dials { namespace model {

  using dials::model::Valid;

  /**
   * A class to hold stuff for an image volume
   */
  class ImageVolume {
  public:

    /**
     * Initialise the class
     * @param frame0 The first frame
     * @param frame1 The last frame
     * @param height The image height
     * @param width The image width
     */
    ImageVolume(int frame0,
                int frame1,
                std::size_t height,
                std::size_t width)
      : frame0_(frame0),
        frame1_(frame1),
        grid_(init_grid(frame0, frame1, height, width)),
        data_(grid_, 0),
        background_(grid_, 0),
        mask_(grid_, 0) {}

    /**
     * Check the arrays all make sense
     */
    bool is_consistent() const {
      return (
        data_.accessor().all_eq(grid_) &&
        background_.accessor().all_eq(grid_) &&
        mask_.accessor().all_eq(grid_));
    }

    /**
     * @returns the first frame
     */
    int frame0() const {
      return frame0_;
    }

    /**
     * @returns The last frame
     */
    int frame1() const {
      return frame1_;
    }

    /**
     * @returns The accessor
     */
    af::c_grid<3> accessor() const {
      return grid_;
    }

    /**
     * @returns The data array
     */
    af::versa < double, af::c_grid<3> > data() const {
      return data_;
    }

    /**
     * @returns The background array
     */
    af::versa < double, af::c_grid<3> > background() const {
      return background_;
    }

    /**
     * @returns The mask array
     */
    af::versa < int, af::c_grid<3> > mask() const {
      return mask_;
    }

    /**
     * Set the image data
     * @param frame The frame number
     * @param data The data array
     * @param mask The mask array
     */
    template <typename T>
    void set_image(
        int frame,
        const af::const_ref< T, af::c_grid<2> > &data,
        const af::const_ref< bool, af::c_grid<2> > &mask) {
      DIALS_ASSERT(frame >= frame0_);
      DIALS_ASSERT(frame < frame1_);
      DIALS_ASSERT(data.accessor().all_eq(mask.accessor()));
      DIALS_ASSERT(data.accessor().all_eq(af::c_grid<2>(grid_[1], grid_[2])));
      std::size_t k = frame - frame0_;
      for (std::size_t j = 0; j < data.accessor()[0]; ++j) {
        for (std::size_t i = 0; i < data.accessor()[1]; ++i) {
          data_(k,j,i) = (double)data(j,i);
          mask_(k,j,i) = (mask(j,i) ? Valid : 0);
        }
      }
    }


  private:

    af::c_grid<3> init_grid(
        int frame0,
        int frame1,
        std::size_t height,
        std::size_t width) const {
      DIALS_ASSERT(frame1 > frame0);
      DIALS_ASSERT(height > 0);
      DIALS_ASSERT(width > 0);
      return af::c_grid<3>(frame1 - frame0, height, width);
    }

    int frame0_;
    int frame1_;
    af::c_grid<3> grid_;
    af::versa <double, af::c_grid<3> > data_;
    af::versa <double, af::c_grid<3> > background_;
    af::versa <int   , af::c_grid<3> > mask_;
  };


  /**
   * A class to hold multiple panel image volumes
   */
  class MultiPanelImageVolume {
  public:

    MultiPanelImageVolume() {}

    /**
     * Add an image volume
     * @param x The image volume
     */
    void add(const ImageVolume &x) {
      if (size() > 0) {
        DIALS_ASSERT(x.frame0() == volume_[0].frame0());
        DIALS_ASSERT(x.frame1() == volume_[0].frame1());
      }
      volume_.push_back(x);
    }

    /**
     * @returns The first frame
     */
    int frame0() const {
      DIALS_ASSERT(size() > 0);
      return volume_[0].frame0();
    }

    /**
     * @returns The last frame
     */
    int frame1() const {
      DIALS_ASSERT(size() > 0);
      return volume_[0].frame1();
    }

    /**
     * Get the image volume for the panel
     * @param index The panel index
     * @returns The image volume
     */
    ImageVolume get(std::size_t index) const {
      DIALS_ASSERT(index < volume_.size());
      return volume_[index];
    }

    /**
     * @returns the number of panels
     */
    std::size_t size() const {
      return volume_.size();
    }

    /**
     * Set the image
     * @param frame
     * @param image The multipanel image
     */
    template <typename T>
    void set_image(int frame, const Image<T> &image) {
      DIALS_ASSERT(image.npanels() == volume_.size());
      for (std::size_t i = 0; i < image.npanels(); ++i) {
        volume_[i].set_image(frame, image.data(i), image.mask(i));
      }
    }

  private:

    af::shared <ImageVolume> volume_;
  };

}};

#endif /* DIALS_MODEL_DATA_IMAGE_VOLUME_H */
