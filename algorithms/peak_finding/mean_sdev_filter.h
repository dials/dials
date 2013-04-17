/*
 * mean_sdev_filter.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_PEAK_FINDING_MEAN_SDEV_FILTER_H
#define DIALS_ALGORITHMS_PEAK_FINDING_MEAN_SDEV_FILTER_H

#include <algorithm>
#include <cmath>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/tiny_types.h>
#include <scitbx/array_family/ref_reductions.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  using scitbx::af::int2;
  using scitbx::af::flex_double;

  /**
   * Calculate the summed area table from the image.
   * @param image The image array
   * @returns The summed area table
   */
  flex_double summed_area_table(const flex_double &image) {

    // Check the dimensions of the image
    DIALS_ASSERT(image.accessor().all().size() == 2);
    
    // Allocate the table
    flex_double table(image.accessor());
    
    // Get the size of the image
    std::size_t ysize = image.accessor().all()[0];
    std::size_t xsize = image.accessor().all()[1];

    // Create the summed area table
    for (std::size_t j = 0; j < ysize; ++j) {
      for (std::size_t i = 0; i < xsize; ++i) {
        double I10 = j > 0 ? table(j - 1, i) : 0;
        double I01 = i > 0 ? table(j, i - 1) : 0;
        double I11 = j > 0 && i > 0 ? table(j - 1, i - 1) : 0; 
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
  flex_double summed_area(const flex_double &image, int2 size) {
    // Check the sizes are valid
    DIALS_ASSERT(image.accessor().all().size() == 2);
    DIALS_ASSERT(size[0] > 0 && size[1] > 0);
    
    // Calculate the summed area table
    flex_double I = summed_area_table(image);

    // Allocate the filtered image
    flex_double sum(image.accessor());
    
    // Get the size of the image
    std::size_t ysize = image.accessor().all()[0];
    std::size_t xsize = image.accessor().all()[1];
    
    // Calculate the local mean at every point
    for (std::size_t j = 0; j < ysize; ++j) {
      for (std::size_t i = 0 ; i < xsize; ++i) {
        int i0 = i - size[0] - 1, i1 = i + size[0];
        int j0 = j - size[1] - 1, j1 = j + size[1];
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

  /**
   * Apply a mean filter to the image
   * @param image The image to filter
   * @param size The size of the filter box (2 * size + 1)
   * @returns The mean filtered image
   */
  flex_double mean_filter(const flex_double &image, int2 size) {
    // Check the sizes are valid
    DIALS_ASSERT(image.accessor().all().size() == 2);
    DIALS_ASSERT(size[0] > 0 && size[1] > 0);
    
    // Calculate the summed area and divide by the number of elements
    flex_double mean = summed_area(image, size);
    double num = (2 * size[0] + 1) * (2 * size[1] + 1);
    for (std::size_t i = 0; i < mean.size(); ++i) {
      mean[i] = mean[i] / num;
    }
    
    // Return the mean filtered image
    return mean;
  }

  /**
   * Apply a sdev filter to the image
   * @param image The image to filter
   * @param mean The mean filtered image
   * @param size The size of the filter box (2 * size + 1)
   * @returns The sdev filtered image
   */
  flex_double sdev_filter(const flex_double &image, const flex_double &mean, 
      int2 size) {
    // Check the sizes are valid
    DIALS_ASSERT(image.accessor().all().size() == 2);
    DIALS_ASSERT(image.accessor().all().all_eq(mean.accessor().all()));
    DIALS_ASSERT(size[0] > 0 && size[1] > 0);

    // Allocate memory for sqr difference image
    flex_double image_sq(image.accessor());
    for (std::size_t i = 0; i < image_sq.size(); ++i) {
      image_sq[i] = image[i] * image[i];
    }
    
    // Calculate and return the sdev filtered image
    flex_double sdev = summed_area(image_sq, size);
    double num = (2 * size[0] + 1) * (2 * size[1] + 1);
    for (std::size_t i = 0; i < sdev.size(); ++i) {
      sdev[i] = sqrt((sdev[i] - mean[i]*mean[i]*num) / num);
    }
    
    // Return the filtered image
    return sdev;
  }
  
  /**
   * Apply a mean and sdev filter to the image and return both images
   * @param image The image filter
   * @param size The window size
   * @returns A pair of mean and sdev filtered images
   */
  std::pair<flex_double, flex_double>
  mean_sdev_filter(const flex_double &image, int2 size) {
    flex_double mean = mean_filter(image, size);
    flex_double sdev = sdev_filter(image, mean, size);
    return std::make_pair(mean, sdev);
  }

}} // namespace dials::algorithms

#endif /* DIALS_ALGORITHMS_PEAK_FINDING_MEAN_SDEV_FILTER_H */

