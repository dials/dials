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
#include <limits>
#include <vector>
#include <boost/unordered_map.hpp>
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

  /**
   * A disjoint-set forest (union-find) for connected component labelling. It
   * uses 4 bytes per element, in contrast to boost::adjacency_list, which
   * needs hundreds of bytes per element for the same job.
   */
  class DisjointSets {
  public:
    DisjointSets() {}

    /** Start with n elements, each in its own set */
    explicit DisjointSets(std::size_t n) : parent_(n) {
      DIALS_ASSERT(n <= std::numeric_limits<int>::max());
      for (std::size_t i = 0; i < n; ++i) {
        parent_[i] = i;
      }
    }

    /** @returns The number of elements */
    std::size_t size() const {
      return parent_.size();
    }

    /**
     * Add an element in a set of its own
     * @returns The index of the new element
     */
    std::size_t add() {
      DIALS_ASSERT(parent_.size() < std::numeric_limits<int>::max());
      parent_.push_back(parent_.size());
      return parent_.size() - 1;
    }

    /**
     * Join the sets containing a and b. The lower index is kept as the root.
     */
    void join(std::size_t a, std::size_t b) {
      std::size_t ra = find(a);
      std::size_t rb = find(b);
      if (ra < rb) {
        parent_[rb] = ra;
      } else if (rb < ra) {
        parent_[ra] = rb;
      }
    }

    /**
     * @returns The label of each element. Sets are numbered in order of their
     * lowest element, which matches the numbering given by
     * boost::connected_components.
     */
    af::shared<int> labels() const {
      af::shared<int> labels(parent_.size(), af::init_functor_null<int>());
      int num = 0;
      for (std::size_t i = 0; i < labels.size(); ++i) {
        // join keeps the lower index as the root, so a parent is never after
        // the element and has already been labelled
        DIALS_ASSERT(parent_[i] <= i);
        labels[i] = (parent_[i] == i) ? num++ : labels[parent_[i]];
      }
      return labels;
    }

  private:
    /** Find the root of a set, halving the path as we go */
    std::size_t find(std::size_t i) {
      while (static_cast<std::size_t>(parent_[i]) != i) {
        parent_[i] = parent_[parent_[i]];
        i = parent_[i];
      }
      return i;
    }

    std::vector<int> parent_;
  };

  template <std::size_t DIM>
  class LabelImageStack;

  /**
   * A class to do connected component labelling on a stack of images.
   */
  template <>
  class LabelImageStack<2> {
  public:
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
    void add_image(const af::const_ref<int, af::c_grid<2> >& image,
                   const af::const_ref<bool, af::c_grid<2> >& mask) {
      // Check the input
      DIALS_ASSERT(image.accessor().all_eq(mask.accessor()));
      DIALS_ASSERT(image.accessor().all_eq(size_));

      // Loop through all the pixels and assign the edges
      std::size_t vertex_a = 0;
      for (std::size_t j = 0; j < size_[0]; ++j) {
        for (std::size_t i = 0; i < size_[1]; ++i) {
          if (mask(j, i)) {
            // Add the vertex
            vertex_a = sets_.add();
            coords_.push_back(vec3<int>(k_, j, i));
            values_.push_back(image(j, i));

            // Add edges to this vertex
            if (i > 0 && mask(j, i - 1)) {
              sets_.join(vertex_a, vertex_a - 1);
            }
            if (j > 0 && mask(j - 1, i)) {
              std::size_t vertex_b = buffer_[i];
              sets_.join(vertex_a, vertex_b - 1);
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
      return sets_.labels();
    }

  private:
    DisjointSets sets_;
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
    void add_image(const af::const_ref<int, af::c_grid<2> >& image,
                   const af::const_ref<bool, af::c_grid<2> >& mask) {
      // Check the input
      DIALS_ASSERT(image.accessor().all_eq(mask.accessor()));
      DIALS_ASSERT(image.accessor().all_eq(size_));

      // Loop through all the pixels and assign the edges
      std::size_t vertex_a = 0;
      for (std::size_t j = 0; j < size_[0]; ++j) {
        for (std::size_t i = 0; i < size_[1]; ++i) {
          if (mask(j, i)) {
            // Add the vertex
            vertex_a = sets_.add();
            coords_.push_back(vec3<int>(k_, j, i));
            values_.push_back(image(j, i));

            // Add edges to this vertex
            if (i > 0 && mask(j, i - 1)) {
              sets_.join(vertex_a, vertex_a - 1);
            }
            if (j > 0 && mask(j - 1, i)) {
              std::size_t vertex_b = buffer_(j - 1, i);
              sets_.join(vertex_a, vertex_b - 1);
            }
            if (k_ > 0 && buffer_(j, i)) {
              std::size_t vertex_b = buffer_(j, i);
              sets_.join(vertex_a, vertex_b - 1);
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
      return sets_.labels();
    }

  private:
    DisjointSets sets_;
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
    void add_pixels(const af::const_ref<int>& values,
                    const af::const_ref<vec3<int> >& coords) {
      DIALS_ASSERT(values.size() == coords.size());
      for (std::size_t i = 0; i < coords.size(); ++i) {
        const vec3<int>& xyz = coords[i];
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
      // Create a hash table of the points
      DisjointSets sets(coords_.size());
      typedef boost::unordered_map<vec3<int>, int, Vec3IntHash> Grid;
      Grid grid(coords_.size(), Vec3IntHash(size_));
      for (std::size_t i = 0; i < coords_.size(); ++i) {
        DIALS_ASSERT(grid.find(coords_[i]) == grid.end());
        grid[coords_[i]] = i;
      }

      // For each point check the pixels to the left in all three dimensions
      // and if they are in the list of pixels then join them. A lookup with
      // find, rather than operator[], avoids inserting the missing neighbours.
      for (std::size_t i = 0; i < coords_.size(); ++i) {
        const vec3<int>& c = coords_[i];
        const vec3<int> left[3] = {vec3<int>(c[0] - 1, c[1], c[2]),
                                   vec3<int>(c[0], c[1] - 1, c[2]),
                                   vec3<int>(c[0], c[1], c[2] - 1)};
        for (std::size_t d = 0; d < 3; ++d) {
          Grid::const_iterator it = grid.find(left[d]);
          if (it != grid.end()) sets.join(i, it->second);
        }
      }

      return sets.labels();
    }

  private:
    // The hash function for the vec3<int> points
    struct Vec3IntHash {
      vec3<int> size_;
      Vec3IntHash(vec3<int> size) : size_(size) {}
      std::size_t operator()(vec3<int> const& x) const {
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
