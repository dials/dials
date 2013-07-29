/*
 * find_overlapping.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_INTEGRATION_FIND_OVERLAPPING_H
#define DIALS_ALGORITHMS_INTEGRATION_FIND_OVERLAPPING_H

#include <vector>
#include <boost/shared_ptr.hpp>
#include <dials/model/data/reflection.h>
#include <dials/model/data/adjacency_list.h>
#include <dials/algorithms/spatial_indexing/detect_collisions.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  using scitbx::af::int6;

  // Helper functions needed for 3D collision detection
  template <> struct bound_coord_type<int6> { typedef int type; };
  template <> int get_minimum_bound<0, int6>(const int6 &b) { return b[0]; }
  template <> int get_maximum_bound<0, int6>(const int6 &b) { return b[1]; }
  template <> int get_minimum_bound<1, int6>(const int6 &b) { return b[2]; }
  template <> int get_maximum_bound<1, int6>(const int6 &b) { return b[3]; }
  template <> int get_minimum_bound<2, int6>(const int6 &b) { return b[4]; }
  template <> int get_maximum_bound<2, int6>(const int6 &b) { return b[5]; }

}} // namespace dials::algorithms

namespace dials { namespace algorithms { namespace shoebox {

  using scitbx::af::int6;
  using dials::model::Reflection;
  using dials::model::ReflectionList;
  using dials::model::AdjacencyList;

  /**
   * Given a set of reflections, find the bounding_boxes that overlap.
   * This function uses a single shot collision detection algorithm to
   * find the colliding bounding_boxes and then puts all the pairs of colliding
   * indices into an adjacency list. Vertices are referred to in the
   * adjacency list by index.
   * @param reflections The reflection list.
   * @returns An adjacency list
   */
  inline
  boost::shared_ptr<AdjacencyList> find_overlapping(
      const ReflectionList &reflections) {

    // Ensure we have a valid number of reflections
    DIALS_ASSERT(reflections.size() > 0);

    // Copy all the reflection bounding_boxes into their own array
    std::vector<int6> bounding_boxes(reflections.size());
    for (std::size_t i = 0; i < reflections.size(); ++i) {
      bounding_boxes[i] = reflections[i].get_bounding_box();
    }

    // Create a list of all the pairs of collisions between bouding boxes.
    std::vector<std::pair<int, int> > collisions;
    detect_collisions3d(bounding_boxes.begin(), bounding_boxes.end(),
                        collisions);

    // Put all the collisions into an adjacency list
    boost::shared_ptr<AdjacencyList> list(new AdjacencyList);
    for (std::size_t i = 0; i < bounding_boxes.size(); ++i) {
      add_vertex(*list);
    }
    for (std::size_t i = 0; i < collisions.size(); ++i) {
      add_edge(collisions[i].first, collisions[i].second, *list);
    }

    // Return the adjacency list
    return list;
  }

}}} // namespace dials::algorithms::shoebox

#endif // DIALS_ALGORITHMS_INTEGRATION_FIND_OVERLAPPING_H
