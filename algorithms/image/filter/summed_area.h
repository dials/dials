/*
 * summed_area.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_IMAGE_FILTER_SUMMED_AREA_H
#define DIALS_ALGORITHMS_IMAGE_FILTER_SUMMED_AREA_H

#include <omptbx/omp_or_stubs.h>

#include <algorithm>
#include <cmath>
#include <scitbx/array_family/tiny_types.h>
#include <scitbx/array_family/ref_reductions.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  using scitbx::af::int2;

  /**
   * Calculate the summed area table from the image.
   * @param image The image array
   * @returns The summed area table
   */
  template <typename T>
  af::versa< T, af::c_grid<2> > summed_area_table(
      const af::const_ref< T, af::c_grid<2> > &image) {

    // Allocate the table
    af::versa< T, af::c_grid<2> > table(image.accessor());

    // Get the size of the image
    std::size_t ysize = image.accessor()[0];
    std::size_t xsize = image.accessor()[1];

    // Create the summed area table
    for (std::size_t j = 0; j < ysize; ++j) {
      for (std::size_t i = 0; i < xsize; ++i) {
        T I10 = j > 0 ? table(j - 1, i) : 0;
        T I01 = i > 0 ? table(j, i - 1) : 0;
        T I11 = j > 0 && i > 0 ? table(j - 1, i - 1) : 0;
        table(j, i) = image(j, i) + I10 + I01 - I11;
      }
    }

    // Return the summed area table
    return table;
  }

  /**
   * Calculate the summed area under each point of the image
   * @param image The image array
   * @param size The size of the rectangle (2 * size + 1)
   * @returns The summed area
   */
  template <typename T>
  af::versa< T, af::c_grid<2> > summed_area(
      const af::const_ref< T, af::c_grid<2> > &image, int2 size) {
    // Check the sizes are valid
    DIALS_ASSERT(size.all_gt(0));

    // Calculate the summed area table
    af::versa< T, af::c_grid<2> > I = summed_area_table<T>(image);

    // Allocate the filtered image
    af::versa< T, af::c_grid<2> > sum(image.accessor());

    // Get the size of the image
    std::size_t ysize = image.accessor()[0];
    std::size_t xsize = image.accessor()[1];

    // Calculate the local mean at every point
    #pragma omp parallel for
    for (std::size_t j = 0; j < ysize; ++j) {
      for (std::size_t i = 0 ; i < xsize; ++i) {
        int i0 = i - size[1] - 1, i1 = i + size[1];
        int j0 = j - size[0] - 1, j1 = j + size[0];
        i1 = i1 < xsize ? i1 : xsize - 1;
        j1 = j1 < ysize ? j1 : ysize - 1;

        double I00 = 0, I10 = 0, I01 = 0, I11 = 0;
        if (i0 >= 0 && j0 >= 0) {
          I00 = I(j0, i0);
          I10 = I(j1, i0);
          I01 = I(j0, i1);
        } else if (i0 >= 0) {
          I10 = I(j1, i0);
        } else if (j0 >= 0) {
          I01 = I(j0, i1);
        }
        I11 = I(j1, i1);

        sum(j, i) = (I11 + I00 - I01 - I10);
      }
    }

    // Return the summed area image
    return sum;
  }

}} // namespace dials::algorithms

#endif /* DIALS_ALGORITHMS_IMAGE_FILTER_SUMMED_AREA_H */
