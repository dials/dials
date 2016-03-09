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
#include <scitbx/array_family/tiny_types.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/model/data/image.h>
#include <dials/model/data/mask_code.h>
#include <dials/error.h>

namespace dials { namespace model {

  using scitbx::af::int6;
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
     * Trim bbox to fit
     */
    int6 trim_bbox(int6 bbox) const {
      int x0 = std::max(0, bbox[0]);
      int y0 = std::max(0, bbox[2]);
      int z0 = std::max(frame0_, bbox[4]);
      int x1 = std::min((int)grid_[2], bbox[1]);
      int y1 = std::min((int)grid_[1], bbox[3]);
      int z1 = std::min(frame1_, bbox[5]);
      DIALS_ASSERT(z1 > z0);
      DIALS_ASSERT(y1 > y0);
      DIALS_ASSERT(x1 > x0);
      return int6(x0, x1, y0, y1, z0, z1);
    }

    /**
     * Extract data with the given bbox
     */
    af::versa < double, af::c_grid<3> > extract_data(int6 bbox) const {
      DIALS_ASSERT(bbox[0] >= 0);
      DIALS_ASSERT(bbox[2] >= 0);
      DIALS_ASSERT(bbox[4] >= frame0_);
      DIALS_ASSERT(bbox[1] <= grid_[2]);
      DIALS_ASSERT(bbox[3] <= grid_[1]);
      DIALS_ASSERT(bbox[5] <= frame1_);
      DIALS_ASSERT(bbox[1] > bbox[0]);
      DIALS_ASSERT(bbox[3] > bbox[2]);
      DIALS_ASSERT(bbox[5] > bbox[4]);
      std::size_t xsize = bbox[1] - bbox[0];
      std::size_t ysize = bbox[3] - bbox[2];
      std::size_t zsize = bbox[5] - bbox[4];
      af::versa< double, af::c_grid<3> > result(
          af::c_grid<3>(zsize, ysize, xsize));
      std::size_t i0 = bbox[0];
      std::size_t j0 = bbox[2];
      std::size_t k0 = bbox[4] - frame0_;
      for (std::size_t k = 0; k < zsize; ++k) {
        for (std::size_t j = 0; j < ysize; ++j) {
          for (std::size_t i = 0; i < xsize; ++i) {
            result(k,j,i) = data_(k+k0, j+j0, i+i0);
          }
        }
      }
      return result;
    }

    /**
     * Extract data with the given bbox
     */
    af::versa < double, af::c_grid<3> > extract_background(int6 bbox) const {
      DIALS_ASSERT(bbox[0] >= 0);
      DIALS_ASSERT(bbox[2] >= 0);
      DIALS_ASSERT(bbox[4] >= frame0_);
      DIALS_ASSERT(bbox[1] <= grid_[2]);
      DIALS_ASSERT(bbox[3] <= grid_[1]);
      DIALS_ASSERT(bbox[5] <= frame1_);
      DIALS_ASSERT(bbox[1] > bbox[0]);
      DIALS_ASSERT(bbox[3] > bbox[2]);
      DIALS_ASSERT(bbox[5] > bbox[4]);
      std::size_t xsize = bbox[1] - bbox[0];
      std::size_t ysize = bbox[3] - bbox[2];
      std::size_t zsize = bbox[5] - bbox[4];
      af::versa< double, af::c_grid<3> > result(
          af::c_grid<3>(zsize, ysize, xsize));
      std::size_t i0 = bbox[0];
      std::size_t j0 = bbox[2];
      std::size_t k0 = bbox[4] - frame0_;
      for (std::size_t k = 0; k < zsize; ++k) {
        for (std::size_t j = 0; j < ysize; ++j) {
          for (std::size_t i = 0; i < xsize; ++i) {
            result(k,j,i) = background_(k+k0, j+j0, i+i0);
          }
        }
      }
      return result;
    }

    /**
     * Extract data with the given bbox
     */
    af::versa < int, af::c_grid<3> > extract_mask(int6 bbox) const {
      DIALS_ASSERT(bbox[0] >= 0);
      DIALS_ASSERT(bbox[2] >= 0);
      DIALS_ASSERT(bbox[4] >= frame0_);
      DIALS_ASSERT(bbox[1] <= grid_[2]);
      DIALS_ASSERT(bbox[3] <= grid_[1]);
      DIALS_ASSERT(bbox[5] <= frame1_);
      DIALS_ASSERT(bbox[1] > bbox[0]);
      DIALS_ASSERT(bbox[3] > bbox[2]);
      DIALS_ASSERT(bbox[5] > bbox[4]);
      std::size_t xsize = bbox[1] - bbox[0];
      std::size_t ysize = bbox[3] - bbox[2];
      std::size_t zsize = bbox[5] - bbox[4];
      af::versa< int, af::c_grid<3> > result(
          af::c_grid<3>(zsize, ysize, xsize));
      std::size_t i0 = bbox[0];
      std::size_t j0 = bbox[2];
      std::size_t k0 = bbox[4] - frame0_;
      for (std::size_t k = 0; k < zsize; ++k) {
        for (std::size_t j = 0; j < ysize; ++j) {
          for (std::size_t i = 0; i < xsize; ++i) {
            result(k,j,i) = mask_(k+k0, j+j0, i+i0);
          }
        }
      }
      return result;
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
