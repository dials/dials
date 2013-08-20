/*
 * connected_components.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_IMAGE_CONNECTED_COMPONENTS_CONNECTED_COMPONENTS_H
#define DIALS_ALGORITHMS_IMAGE_CONNECTED_COMPONENTS_CONNECTED_COMPONENTS_H

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <scitbx/vec3.h>
#include <scitbx/array_family/tiny_types.h>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/shared.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  using scitbx::vec3;
  using scitbx::af::int2;
  using scitbx::af::flex_int;
  using scitbx::af::flex_bool;
  using scitbx::af::flex_double;
  using scitbx::af::flex_grid;

  /**
   * A class to do connected component labelling on a stack of images.
   */
  class LabelImageStack {
  public:

    // Adjacency list type
    typedef boost::adjacency_list<
      boost::listS,
      boost::vecS,
      boost::undirectedS> AdjacencyList;

    /**
     * Initialise the class with the size of the desired image.
     * @param size The size of the images
     */
    LabelImageStack(int2 size)
      : buffer_(flex_grid<>(size[0], size[1])),
        size_(size),
        k_(0) {}

    /**
     * Add another image to be labelled
     * @param image The image to use
     * @param mask The mask to use
     */
    void add_image(const flex_double &image, const flex_bool &mask) {

      // Check the input
      DIALS_ASSERT(mask.accessor().all().size() == 2);
      DIALS_ASSERT(mask.accessor().all()[0] == size_[0]);
      DIALS_ASSERT(mask.accessor().all()[1] == size_[1]);

      // Loop through all the pixels and assign the edges
      std::size_t vertex_a = 0;
      for (std::size_t j = 0; j < size_[0]; ++j) {
        for (std::size_t i = 0; i < size_[1]; ++i) {
          if (mask(j, i)) {

            // Add the vertex
            vertex_a = add_vertex(graph_);
            coords_.push_back(vec3<int>(k_, j, i));
            values_.push_back(image(j, i));

            // Add edges to this vertex
            if (i > 0 && mask(j, i - 1)) {
              boost::add_edge(vertex_a, vertex_a - 1, graph_);
            }
            if (j > 0 && mask(j - 1, i)) {
              std::size_t vertex_b = buffer_(j - 1, i);
              boost::add_edge(vertex_a, vertex_b - 1, graph_);
            }
            if (k_ > 0 && buffer_(j, i)) {
              std::size_t vertex_b = buffer_(j, i);
              boost::add_edge(vertex_a, vertex_b - 1, graph_);
            }
            buffer_(j, i) = vertex_a + 1;
          } else {
            buffer_(j, i) = 0;
          }
        }
      }

      // Increment image number
      k_++;
    }

    /**
     * @returns The list of valid point coordinates
     */
    scitbx::af::shared< vec3<int> > coords() const {
      return coords_;
    }

    /**
     * @returns The list of valid point values
     */
    scitbx::af::shared<double> values() const {
      return values_;
    }

    /**
     * Do the connected component labelling and get the labels for each
     * of the good pixels given to the algorithm.
     * @returns The list of labels
     */
    scitbx::af::shared<int> labels() const {
      scitbx::af::shared<int> labels(num_vertices(graph_));
      int num = boost::connected_components(graph_, &labels[0]);
      DIALS_ASSERT(num <= labels.size());
      return labels;
    }

  private:

    AdjacencyList graph_;
    scitbx::af::shared< vec3<int> > coords_;
    scitbx::af::shared<double> values_;
    scitbx::af::flex_size_t buffer_;
    int2 size_;
    std::size_t k_;
  };

}} // namespace dials::algorithms

#endif /* DIALS_ALGORITHMS_IMAGE_CONNECTED_COMPONENTS_CONNECTED_COMPONENTS_H */
