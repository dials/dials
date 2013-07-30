/*
 * masker.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_INTEGRATION_MASKER_H
#define DIALS_ALGORITHMS_INTEGRATION_MASKER_H

#include <scitbx/vec3.h>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/tiny_types.h>
#include <dials/model/data/reflection.h>
#include <dials/model/data/adjacency_list.h>
#include <dials/algorithms/shoebox/mask_code.h>
#include <dials/error.h>

namespace dials { namespace algorithms { namespace shoebox {

  using scitbx::vec3;
  using scitbx::af::int6;
  using scitbx::af::flex_bool;
  using scitbx::af::flex_int;
  using dials::model::Reflection;
  using dials::model::ReflectionList;
  using dials::model::AdjacencyList;

  /** Class to calculate the shoebox masks for all reflections */
  class Masker {
  public:

    // Useful typedefs
    typedef AdjacencyList::adjacency_iterator adjacency_iterator;
    typedef std::pair<adjacency_iterator,adjacency_iterator>
      adjacency_iterator_range;

    /**
     * Initialise the algorithm with the detector mask
     * @param detector_mask The mask of good/bad detector pixels.
     */
    Masker(const flex_bool &detector_mask)
      : detector_mask_(detector_mask),
        detector_size_(detector_mask.accessor().all()) {
      DIALS_ASSERT(detector_size_.size() == 2);
      DIALS_ASSERT(detector_size_[0] > 0 && detector_size_[1] > 0);
    }

    /**
     * The entry point of the functor. The list of reflections is queried
     * against the adjacency list to find overlapping reflections. Those
     * that have no overlaps their whole mask set to 1 to indicate that
     * the reflection owns all it's pixels. Those that do overlap are
     * then compared pixel-by-pixel in the overlapping region and the
     * reflection whose predicted central location is closer to the pixel
     * will gain ownership of the pixel (mask value 1).
     *
     * @param reflections The list of reflections
     * @param adjacency_list The adjacency_list
     */
    void operator()(ReflectionList &reflections,
        const boost::shared_ptr<AdjacencyList> &adjacency_list) const {

      // Begin by setting all the mask values to 1
      initialise_mask(reflections);

      // Loop through all the reflections
      for (std::size_t i = 0; i < reflections.size(); ++i) {

        // Get a reference to the reflection
        Reflection &r = reflections[i];

        // Get the list of overlapping reflectionss
        adjacency_iterator_range range = adjacent_vertices(i, *adjacency_list);
        for (adjacency_iterator it = range.first; it != range.second; ++it) {
          if (i < *it) {
            assign_ownership(r, reflections[*it]);
          }
        }
      }
    }

  private:

    /**
     * The distance between two points
     * @param a Point a
     * @param b Point b
     * @returns The distance from a to b
     */
    double distance(vec3<double> a, vec3<double> b) const {
      return (b - a).length_sq();
    }

    /**
     * Get a coordinate from the reflection.
     * @param r The reflection
     * @returns An (x, y, z) coordinate
     */
    vec3<double> reflection_coord(Reflection &r) const {
      return vec3<double>(
        r.get_image_coord_px()[0],
        r.get_image_coord_px()[1],
        r.get_frame_number());
    }

    /**
     * Get a coordinate from the voxel.
     * @param i The voxel x index
     * @param j The voxel y index
     * @param k THe voxel z index
     * @returns An (x, y, z) coordinate
     */
    vec3<double> voxel_coord(double i, double j, double k) const {
      return vec3<double>(i+0.5, j+0.5, k+0.5);
    }

    /**
     * Assign ownership of all the pixels in the overlapping range.
     * @param a Reflection a
     * @param b Reflection b
     * @throws RuntimeError if reflections to do overlap.
     */
    void assign_ownership(Reflection &a, Reflection &b) const {

      // Get the reflection mask arrays
      flex_int mask_a = a.get_shoebox_mask();
      flex_int mask_b = b.get_shoebox_mask();

      // Get the sizes of the masks
      flex_int::index_type size_a = mask_a.accessor().all();
      flex_int::index_type size_b = mask_b.accessor().all();

      // Get the bounding boxes
      int6 bbox_a = a.get_bounding_box();
      int6 bbox_b = b.get_bounding_box();

      // Get the reflection coordinates in pixels
      vec3<double> coord_a = reflection_coord(a);
      vec3<double> coord_b = reflection_coord(b);

      // Get the status for each reflection
      bool a_status = a.is_valid();
      bool b_status = b.is_valid();

      // Get range to iterate over
      int i0 = std::max(bbox_a[0], bbox_b[0]);
      int i1 = std::min(bbox_a[1], bbox_b[1]);
      int j0 = std::max(bbox_a[2], bbox_b[2]);
      int j1 = std::min(bbox_a[3], bbox_b[3]);
      int k0 = std::max(bbox_a[4], bbox_b[4]);
      int k1 = std::min(bbox_a[5], bbox_b[5]);

      // Ensure ranges are valid
      DIALS_ASSERT(k1 > k0 && j1 > j0 && i1 > i0);

      if (a_status == true) {
        DIALS_ASSERT(i0 - bbox_a[0] >= 0 && i1 - bbox_a[0] <= size_a[2]);
        DIALS_ASSERT(j0 - bbox_a[2] >= 0 && j1 - bbox_a[2] <= size_a[1]);
        DIALS_ASSERT(k0 - bbox_a[4] >= 0 && k1 - bbox_a[4] <= size_a[0]);
      }
      if (b_status == true) {
        DIALS_ASSERT(i0 - bbox_b[0] >= 0 && i1 - bbox_b[0] <= size_b[2]);
        DIALS_ASSERT(j0 - bbox_b[2] >= 0 && j1 - bbox_b[2] <= size_b[1]);
        DIALS_ASSERT(k0 - bbox_b[4] >= 0 && k1 - bbox_b[4] <= size_b[0]);
      }

      // Iterate over range of indices
      for (int k = k0; k < k1; ++k) {
        for (int j = j0; j < j1; ++j) {
          for (int i = i0; i < i1; ++i) {

            // Set the coordinate
            vec3<double> coord_c(voxel_coord(i, j, k));

            // Set the mask indices
            int ka = k - bbox_a[4], kb = k - bbox_b[4];
            int ja = j - bbox_a[2], jb = j - bbox_b[2];
            int ia = i - bbox_a[0], ib = i - bbox_b[0];

            // If the distance from a to c is less than b to c then
            // set the mask for a to 1 and b to 0, otherwise set
            // b to 1 and a to 0.
            if (distance(coord_a, coord_c) < distance(coord_b, coord_c)) {
              if (b_status == true) {
                mask_b(kb, jb, ib) = 0;
              }
            } else {
              if (a_status == true) {
                mask_a(ka, ja, ia) = 0;
              }
            }
          }
        }
      }
    }

    /**
     * Set all the shoebox mask values for the reflection to value
     * @param reflection The reflection
     */
    void initialise_mask(Reflection &reflection) const {

      // Get the mask and roi from the reflection
      flex_int mask = reflection.get_shoebox_mask();
      flex_int::index_type size = mask.accessor().all();
      int6 roi = reflection.get_bounding_box();

      // Check sizes are ok
      DIALS_ASSERT(size.size() == 3);
      DIALS_ASSERT(size[2] == (roi[1] - roi[0]));
      DIALS_ASSERT(size[1] == (roi[3] - roi[2]));
      DIALS_ASSERT(size[0] == (roi[5] - roi[4]));

      // Loop through all the pixels in the reflection mask and check them
      // against the detector mask. If the detector mask pixel is ok then
      // set the reflection mask pixel to the given value, otherwise set it
      // to zero.
      for (std::size_t j = 0; j < size[1]; ++j) {
        for (std::size_t i = 0; i < size[2]; ++i) {
          int di = roi[0] + i;
          int dj = roi[2] + j;
          bool mask_value = false;
          if (di >= 0 && di < detector_size_[1] &&
              dj >= 0 && dj < detector_size_[0]) {
            mask_value = detector_mask_(dj, di);
          }
          for (std::size_t k = 0; k < size[0]; ++k) {
            mask(k, j, i) = (mask_value ? Valid : 0);
          }
        }
      }
    }

    /**
     * Set all the shoebox mask values for all the reflections to value
     * @param reflections The reflection list
     */
    void initialise_mask(ReflectionList &reflections) const {
      for (std::size_t i = 0; i < reflections.size(); ++i) {
        if (reflections[i].is_valid()) {
          initialise_mask(reflections[i]);
        }
      }
    }

    flex_bool detector_mask_;
    flex_bool::index_type detector_size_;
  };

}}} // namespace dials::algorithms::shoebox

#endif /* DIALS_ALGORITHMS_INTEGRATION_MASKER_H */
