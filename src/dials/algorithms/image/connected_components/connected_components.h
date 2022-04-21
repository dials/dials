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

#include <ctime>
#include <algorithm>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include <scitbx/array_family/tiny_types.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  using scitbx::vec2;
  using scitbx::vec3;
  using scitbx::af::int2;
  using scitbx::af::int3;

  template <std::size_t DIM>
  class LabelImageStack;

  /**
   * A class to do connected component labelling on a stack of images.
   */
  template <>
  class LabelImageStack<2> {
  public:
    // Adjacency list type
    typedef boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS>
      AdjacencyList;

    /**
     * Initialise the class with the size of the desired image.
     * @param size The size of the images
     */
    LabelImageStack(int2 size)
        : buffer_(size[1], af::init_functor_null<std::size_t>()), size_(size), k_(0) {}

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
    void add_image(const af::const_ref<int, af::c_grid<2> > &image,
                   const af::const_ref<bool, af::c_grid<2> > &mask) {
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
    af::shared<vec3<int> > coords() const {
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
      af::shared<int> labels(num_vertices(graph_), af::init_functor_null<int>());
      int num = boost::connected_components(graph_, &labels[0]);
      DIALS_ASSERT(num <= labels.size());
      return labels;
    }

  private:
    AdjacencyList graph_;
    af::shared<vec3<int> > coords_;
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
    typedef boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS>
      AdjacencyList;

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
    void add_image(const af::const_ref<int, af::c_grid<2> > &image,
                   const af::const_ref<bool, af::c_grid<2> > &mask) {
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
    af::shared<vec3<int> > coords() const {
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
      af::shared<int> labels(num_vertices(graph_), af::init_functor_null<int>());
      int num = boost::connected_components(graph_, &labels[0]);
      DIALS_ASSERT(num <= labels.size());
      return labels;
    }

  private:
    AdjacencyList graph_;
    af::shared<vec3<int> > coords_;
    af::shared<int> values_;
    af::versa<std::size_t, af::c_grid<2> > buffer_;
    int2 size_;
    std::size_t k_;
  };

  /**
   * Class to do connected component labelling of input pixels and coords
   */
  class LabelPixels {
  public:
    /** @param size Size of the 3D volume (z, y, x) */
    LabelPixels(int3 size) : size_(size) {}

    /** @returns The size of the 3D volume */
    int3 size() const {
      return size_;
    }

    /**
     * Add pixels to be labelled
     * @param values The pixel values
     * @param coords The pixel coords
     */
    void add_pixels(const af::const_ref<int> &values,
                    const af::const_ref<vec3<int> > &coords) {
      DIALS_ASSERT(values.size() == coords.size());
      for (std::size_t i = 0; i < coords.size(); ++i) {
        const vec3<int> &xyz = coords[i];
        DIALS_ASSERT(xyz[0] >= 0 && xyz[0] < size_[2]);
        DIALS_ASSERT(xyz[1] >= 0 && xyz[1] < size_[1]);
        DIALS_ASSERT(xyz[2] >= 0 && xyz[2] < size_[0]);
        vec3<int> zyx(xyz[2], xyz[1], xyz[0]);
        coords_.push_back(zyx);
        values_.push_back(values[i]);
      }
    }

    /**
     * @returns The list of valid point coordinates
     */
    af::shared<vec3<int> > coords() const {
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
      // Adjacency list type
      typedef boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS>
        AdjacencyList;

      // Create a hash table of the points
      AdjacencyList graph(coords_.size());
      boost::unordered_map<vec3<int>, int, Vec3IntHash> grid(coords_.size(),
                                                             Vec3IntHash(size_));
      for (std::size_t i = 0; i < coords_.size(); ++i) {
        DIALS_ASSERT(grid[coords_[i]] == 0);
        grid[coords_[i]] = i + 1;
      }

      // For each point check the pixels to the left in all three dimensins
      // and if they are in the list of pixels then add the edges.
      for (std::size_t i = 0; i < coords_.size(); ++i) {
        int j;
        j = grid[vec3<int>(coords_[i][0] - 1, coords_[i][1], coords_[i][2])];
        if (j != 0) boost::add_edge(i, j - 1, graph);
        j = grid[vec3<int>(coords_[i][0], coords_[i][1] - 1, coords_[i][2])];
        if (j != 0) boost::add_edge(i, j - 1, graph);
        j = grid[vec3<int>(coords_[i][0], coords_[i][1], coords_[i][2] - 1)];
        if (j != 0) boost::add_edge(i, j - 1, graph);
      }

      // Do the connected components
      af::shared<int> labels(num_vertices(graph), af::init_functor_null<int>());
      int num = boost::connected_components(graph, &labels[0]);
      DIALS_ASSERT(num <= labels.size());
      return labels;
    }

  private:
    // The hash function for the vec3<int> points
    struct Vec3IntHash {
      vec3<int> size_;
      Vec3IntHash(vec3<int> size) : size_(size) {}
      std::size_t operator()(vec3<int> const &x) const {
        boost::hash<int> hasher;
        return hasher(x[0] + x[1] * size_[2] + x[2] * size_[1] * size_[2]);
      }
    };

    af::shared<vec3<int> > coords_;
    af::shared<int> values_;
    int3 size_;
  };

}}  // namespace dials::algorithms

#endif /* DIALS_ALGORITHMS_IMAGE_CONNECTED_COMPONENTS_CONNECTED_COMPONENTS_H */
