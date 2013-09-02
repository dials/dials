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

#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/tiny_types.h>
#include <scitbx/array_family/small.h>
#include <dials/error.h>

namespace dials { namespace model {

  using scitbx::af::int2;
  using scitbx::af::int6;
  using scitbx::af::small;
  using scitbx::af::flex_double;
  using scitbx::af::flex_int;
  using scitbx::af::flex_bool;
  using scitbx::af::flex_grid;

  /**
   * An enumeration of shoebox mask codes. If a pixel is labelled as:
   *  a) Valid. This means that the pixel belongs to this reflection and
   *     is not a bad pixel etc.
   *  b) Background. This means that the pixel is to be used for background
   *     determination.
   *  c) Foreground. This means that the pixel is to be used for integration.
   *  d) Strong. This means that the pixel is defined as strong
   */
  enum MaskCode {
    Valid = (1 << 0),       ///< Pixel is valid for this shoebox
    Background = (1 << 1),  ///< Pixel is in the background
    Foreground = (1 << 2),  ///< Pixel is in the foreground
    Strong = (1 << 3),      ///< Pixel is a strong pixel
  };

  /**
   * A class to hold shoebox information
   */
  struct Shoebox {

    int6 bbox;          ///< The bounding box
    flex_double data;   ///< The shoebox data
    flex_int mask;      ///< The shoebox mask

    /**
     * Initialise the shoebox
     */
    Shoebox()
      : bbox(0, 0, 0, 0, 0, 0),
        data(flex_grid<>(0, 0, 0)),
        mask(flex_grid<>(0, 0, 0)) {}

    /**
     * Initialise the shoebox
     * @param bbox_ The bounding box to initialise with
     */
    Shoebox(const int6 &bbox_)
      : bbox(bbox_),
        data(flex_grid<>(0, 0, 0)),
        mask(flex_grid<>(0, 0, 0)) {}

    /**
     * Allocate the mask and data from the bounding box
     */
    void allocate() {
      data = flex_double(flex_grid<>(zsize(), ysize(), xsize()));
      mask = flex_int(flex_grid<>(zsize(), ysize(), xsize()));
    }

    /**
     * Deallocate the mask and data arrays
     */
    void deallocate() {
      data = flex_double(flex_grid<>(0, 0, 0));
      mask = flex_int(flex_grid<>(0, 0, 0));
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
      return bbox[1] - bbox[0];
    }

    /** @returns The y size */
    std::size_t ysize() const {
      DIALS_ASSERT(bbox[3] >= bbox[2]);
      return bbox[3] - bbox[2];
    }

    /** @returns The z size */
    std::size_t zsize() const {
      DIALS_ASSERT(bbox[5] >= bbox[4]);
      return bbox[5] - bbox[4];
    }

    /** @returns The offset */
    small<long, 10> offset() const {
      small<long, 10> off(3);
      off[0] = zoffset();
      off[1] = yoffset();
      off[2] = xoffset();
      return off;
    }

    /** @returns The size */
    small<long, 10> size() const {
      small<long, 10> sz(3);
      sz[0] = zsize();
      sz[1] = ysize();
      sz[2] = xsize();
      return sz;
    }

    /** @return True/False whether the array and bbox sizes are consistent */
    bool is_consistent() const {
      bool result = true;
      result = result && (data.accessor().all().size() == 3);
      result = result && (mask.accessor().all().size() == 3);
      result = result && (data.accessor().all().all_eq(size()));
      result = result && (mask.accessor().all().all_eq(size()));
      return result;
    }

    /**
     * Check if the bounding box has points outside the image range.
     * @param image_size The image size
     * @param scan_range The scan range
     * @returns True/False
     */
    bool is_bbox_within_image_volume(small<long,10> image_size,
        int2 scan_range) const {
      DIALS_ASSERT(image_size.size() == 2);
      return bbox[0] >= 0 && bbox[1] < image_size[1] &&
             bbox[2] >= 0 && bbox[3] < image_size[0] &&
             bbox[4] >= scan_range[0] && bbox[5] < scan_range[1];
    }

    /**
     * Check if the bounding box has points that cover bad pixels
     * @param mask The mask array
     * @returns True/False
     */
    inline
    bool does_bbox_contain_bad_pixels(const flex_bool &mask) const {
      DIALS_ASSERT(mask.accessor().all().size() == 2);
      std::size_t ysize = mask.accessor().all()[0];
      std::size_t xsize = mask.accessor().all()[1];
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
        if (mask[i] & code) {
          count++;
        }
      }
      return count;
    }

    /**
     * Test to see if shoeboxs contain the same data
     * @param rhs The other shoebox
     * @returns True/False. They are the same
     */
    bool operator==(const Shoebox &rhs) const {
      return ((bbox.all_eq(rhs.bbox)) &&
              (data.all_eq(rhs.data)) &&
              (mask.all_eq(rhs.mask)));
    }

    /**
     * Test to see if shoeboxs contain the same data
     * @param rhs The other shoebox
     * @returns True/False. They are not the same
     */
    bool operator!=(const Shoebox &rhs) const {
      return !(*this == rhs);
    }
  };

}};

#endif /* DIALS_MODEL_DATA_SHOEBOX_H */
