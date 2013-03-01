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
#ifndef DIALS_ALGORITHMS_INTEGRATION_SHOEBOX_H
#define DIALS_ALGORITHMS_INTEGRATION_SHOEBOX_H

#include <scitbx/array_family/tiny_types.h>

namespace dials { namespace algorithms {

  using scitbx::af::int2;
  using scitbx::af::int6;

  /**
   * A 3D shoebox class which contains data in the form of a 6 element array
   * like: [minx, maxx, miny, maxy, minz, maxz]
   */
  class Shoebox3d : int6 {
  public:

    /** Default initialisation to 0 */
    Shoebox3d() : int6(0, 0, 0, 0, 0, 0) {}

    /** Initialise shoebox from int6 */
    Shoebox3d(int6 a) : int6(a) {}

    /** Initialise shoebox from min/max x/y/z */
    Shoebox3d(int x0, int x1, int y0, int y1, int z0, int z1)
      : int6(x0, x1, y0, y1, z0, z1) {}

    /** Get the min x value */
    int min_x() const {
      return elems[0];
    }

    /** Get the max x value */
    int max_x() const {
      return elems[1];
    }

    /** Get the min y value */
    int min_y() const {
      return elems[2];
    }

    /** Get the max y value */
    int max_y() const {
      return elems[3];
    }

    /** Get the min z value */
    int min_z() const {
      return elems[4];
    }

    /** Get the max z value */
    int max_z() const {
      return elems[5];
    }

    /** Get the range of x values */
    int2 range_x() const {
      return int2(min_x(), max_x());
    }

    /** Get the range of y values */
    int2 range_y() const {
      return int2(min_y(), max_y());
    }

    /** Get the range of z values */
    int2 range_z() const {
      return int2(min_z(), max_z());
    }
  };

}} // namespace dials::algorithms

#endif // DIALS_ALGORITHMS_INTEGRATION_SHOEBOX_H
