/*
 * basic_shoebox.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_MODEL_DATA_BASIC_SHOEBOX_H
#define DIALS_MODEL_DATA_BASIC_SHOEBOX_H

#include <scitbx/array_family/tiny_types.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/error.h>

namespace dials { namespace model {

  using scitbx::af::int2;
  using scitbx::af::int3;
  using scitbx::af::int6;

  class BasicShoebox {
  public:

    std::size_t panel;
    int6 bbox;
    af::versa< int, af::c_grid<3> > data;

    /**
     * Initialise the shoebox
     */
    BasicShoebox()
      : panel(0),
        bbox(0, 0, 0, 0, 0, 0) {}

    /**
     * Initialise the shoebox
     * @param panel_ The panel number
     * @param bbox_ The bounding box to initialise with
     */
    BasicShoebox(std::size_t panel_, const int6 &bbox_)
      : panel(panel_),
        bbox(bbox_) {}

    /**
     * Allocate the mask and data from the bounding box
     */
    void allocate() {
      af::c_grid<3> accessor(zsize(), ysize(), xsize());
      data = af::versa< int, af::c_grid<3> >(accessor, 0);
    }

    /**
     * Deallocate the mask and data arrays
     */
    void deallocate() {
      af::c_grid<3> accessor(0, 0, 0);
      data = af::versa< int, af::c_grid<3> >(accessor);
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
    int3 offset() const {
      return int3(zoffset(), yoffset(), xoffset());
    }

    /** @returns The size */
    int3 size() const {
      return int3(zsize(), ysize(), xsize());
    }

    /** @return True/False whether the array and bbox sizes are consistent */
    bool is_consistent() const {
      return data.accessor().all_eq(size());
    }

    /**
     * Test to see if shoeboxs contain the same data
     * @param rhs The other shoebox
     * @returns True/False. They are the same
     */
    bool operator==(const BasicShoebox &rhs) const {
      return ((bbox.all_eq(rhs.bbox)) &&
              (data.all_eq(rhs.data)));
    }

    /**
     * Test to see if shoeboxs contain the same data
     * @param rhs The other shoebox
     * @returns True/False. They are not the same
     */
    bool operator!=(const BasicShoebox &rhs) const {
      return !(*this == rhs);
    }
  };

}} // namespace dials::model

#endif /* DIALS_MODEL_DATA_BASIC_SHOEBOX_H */
