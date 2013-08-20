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

  using scitbx::af::int6;
  using scitbx::af::small;
  using scitbx::af::flex_double;
  using scitbx::af::flex_int;
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
     * Allocate the mask and data from the bounding box
     */
    void allocate() {
      flex_grid<> accessor(zsize(), ysize(), xsize());
      data.resize(accessor);
      mask.resize(accessor);
    }

    /**
     * Deallocate the mask and data arrays
     */
    void deallocate() {
      flex_grid<> accessor(0, 0, 0);
      data.resize(accessor);
      mask.resize(accessor);
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
    bool consistent() const {
      bool result = true;
      result = result && (data.accessor().all().size() == 3);
      result = result && (mask.accessor().all().size() == 3);
      result = result && (data.accessor().all().all_eq(size()));
      result = result && (mask.accessor().all().all_eq(size()));
      return result;
    }
  };

}};

#endif /* DIALS_MODEL_DATA_SHOEBOX_H */
