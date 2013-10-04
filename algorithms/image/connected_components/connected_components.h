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
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  using scitbx::vec3;
  using scitbx::af::int2;

  template <std::size_t DIM>
  class LabelImageStack;

  /**
   * A class to do connected component labelling on a stack of images.
   */
  template <>
  class LabelImageStack<2> {
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
      : buffer_(size[1], af::init_functor_null<std::size_t>()),
        size_(size),
        k_(0) {}

    /**
     * @returns The image size
     */
    int2 size() const {
      return size_;
    }

    /**
     * @returns The number of images processed
     */
    int num_images() const {
      return k_;
    }

    /**
     * Add another image to be labelled
     * @param image The image to use
     * @param mask The mask to use
     */
    void add_image(const af::const_ref< int, af::c_grid<2> > &image,
                   const af::const_ref< bool, af::c_grid<2> > &mask) {

      // Check the input
      DIALS_ASSERT(image.accessor().all_eq(mask.accessor()));
      DIALS_ASSERT(image.accessor().all_eq(size_));

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
              std::size_t vertex_b = buffer_[i];
              boost::add_edge(vertex_a, vertex_b - 1, graph_);
            }
            buffer_[i] = vertex_a + 1;
          } else {
            buffer_[i] = 0;
          }
        }
      }

      // Increment image number
      k_++;
    }

    /**
     * @returns The list of valid point coordinates
     */
    af::shared< vec3<int> > coords() const {
      return coords_;
    }

    /**
     * @returns The list of valid point values
     */
    af::shared<int> values() const {
      return values_;
    }

    /**
     * Do the connected component labelling and get the labels for each
     * of the good pixels given to the algorithm.
     * @returns The list of labels
     */
    af::shared<int> labels() const {
      af::shared<int> labels(num_vertices(graph_),
        af::init_functor_null<int>());
      int num = boost::connected_components(graph_, &labels[0]);
      DIALS_ASSERT(num <= labels.size());
      return labels;
    }

  private:

    AdjacencyList graph_;
    af::shared< vec3<int> > coords_;
    af::shared<int> values_;
    af::shared<std::size_t> buffer_;
    int2 size_;
    std::size_t k_;
  };


  /**
   * A class to do connected component labelling on a stack of images.
   */
  template <>
  class LabelImageStack<3> {
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
      : buffer_(af::c_grid<2>(size), af::init_functor_null<std::size_t>()),
        size_(size),
        k_(0) {}

    /**
     * @returns The image size
     */
    int2 size() const {
      return size_;
    }

    /**
     * @returns The number of images processed
     */
    int num_images() const {
      return k_;
    }

    /**
     * Add another image to be labelled
     * @param image The image to use
     * @param mask The mask to use
     */
    void add_image(const af::const_ref< int, af::c_grid<2> > &image,
                   const af::const_ref< bool, af::c_grid<2> > &mask) {

      // Check the input
      DIALS_ASSERT(image.accessor().all_eq(mask.accessor()));
      DIALS_ASSERT(image.accessor().all_eq(size_));

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
    af::shared< vec3<int> > coords() const {
      return coords_;
    }

    /**
     * @returns The list of valid point values
     */
    af::shared<int> values() const {
      return values_;
    }

    /**
     * Do the connected component labelling and get the labels for each
     * of the good pixels given to the algorithm.
     * @returns The list of labels
     */
    af::shared<int> labels() const {
      af::shared<int> labels(num_vertices(graph_),
        af::init_functor_null<int>());
      int num = boost::connected_components(graph_, &labels[0]);
      DIALS_ASSERT(num <= labels.size());
      return labels;
    }

  private:

    AdjacencyList graph_;
    af::shared< vec3<int> > coords_;
    af::shared<int> values_;
    af::versa< std::size_t, af::c_grid<2> > buffer_;
    int2 size_;
    std::size_t k_;
  };

}} // namespace dials::algorithms

#endif /* DIALS_ALGORITHMS_IMAGE_CONNECTED_COMPONENTS_CONNECTED_COMPONENTS_H */
