/*
 * label_pixels.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_PEAK_FINDING_LABEL_PIXELS_H
#define DIALS_ALGORITHMS_PEAK_FINDING_LABEL_PIXELS_H

#include <boost/unordered_map.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <scitbx/vec3.h>
#include <dials/array_family/scitbx_shared_and_versa.h>

namespace dials { namespace algorithms {

  using scitbx::vec3;

  /** A class to label the pixels into connected components. */
  class LabelPixels {
  public:

    // Adjacency list type
    typedef boost::adjacency_list<
      boost::listS,
      boost::vecS,
      boost::undirectedS> AdjacencyList;

    /** Construct the class with the grid size */
    LabelPixels(vec3<int> grid_size)
      : grid_size_(grid_size) {}

    /**
     * Calculate the connected components
     * @param pixels The list of pixel coordinates
     * @returns A list of labels into the pixel array.
     */
    af::shared<int> operator()(const af::const_ref< vec3<int> > &pixels) {
      af::shared<int> labels(pixels.size());

      // Build the adjacency list from the pixels
      AdjacencyList graph(pixels.size());
      build_adjacency_list(pixels, graph);

      // Calculate the connected components
      int num = boost::connected_components(graph, &labels[0]);

      // Return the connected components
      return labels;
    }

  private:

    // The hash function for the vec3<int> points
    struct Vec3IntHash {
      vec3<int> size_;
      Vec3IntHash(vec3<int> size) : size_(size) {}
      std::size_t operator()(vec3<int> const &x) const {
        boost::hash<int> hasher;
        return hasher(x[0] + x[1]*size_[0] + x[2]*size_[0]*size_[1]);
      }
    };

    /**
     * Build the adjacency list from the set of points
     * @param pixels The pixel list
     * @param graph The adjacency list
     */
    void build_adjacency_list(const af::const_ref< vec3<int> > &pixels,
        AdjacencyList &graph) {

      // Create a hash table of the points
      boost::unordered_map<vec3<int>, int, Vec3IntHash> grid(
        pixels.size(), Vec3IntHash(grid_size_));
      for (std::size_t i = 0; i < pixels.size(); ++i) {
        grid[pixels[i]] = i + 1;
      }

      // For each point check the pixels to the left in all three dimensins
      // and if they are in the list of pixels then add the edges.
      for (std::size_t i = 0; i < pixels.size(); ++i) {
        int j;
        j = grid[vec3<int>(pixels[i][0]-1, pixels[i][1], pixels[i][2])];
        if (j != 0) {
          add_edge(i, j-1, graph);
        }
        j = grid[vec3<int>(pixels[i][0], pixels[i][1]-1, pixels[i][2])];
        if (j != 0) {
          add_edge(i, j-1, graph);
        }
        j = grid[vec3<int>(pixels[i][0], pixels[i][1], pixels[i][2]-1)];
        if (j != 0) {
          add_edge(i, j-1, graph);
        }
      }
    }

    vec3<int> grid_size_;
  };

}} // namespace dials::algorithms

#endif /* DIALS_ALGORITHMS_PEAK_FINDING_LABEL_PIXELS_H */
