/*
 * pixel_list.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_MODEL_DATA_PIXEL_LIST_H
#define DIALS_MODEL_DATA_PIXEL_LIST_H

#include <scitbx/vec3.h>
#include <scitbx/array_family/tiny_types.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/algorithms/image/connected_components/connected_components.h>
#include <dials/error.h>

namespace dials { namespace model {

  using scitbx::vec3;
  using scitbx::af::int2;

  namespace detail {

    inline bool lessthan(const vec3<int> &a, const vec3<int> &b) {
      return (
        a[0] < b[0]
          ? true
          : (a[0] == b[0]
               ? (a[1] < b[1] ? true
                              : (a[1] == b[1] ? (a[2] < b[2] ? true : false) : false))
               : false));
    }

  }  // namespace detail

  /**
   * A class to hold a list of pixels
   */
  class PixelList {
    int frame_;
    int2 size_;
    af::shared<double> value_;
    af::shared<std::size_t> index_;

  public:
    /**
     * Initialise an empty list
     */
    PixelList() : frame_(0), size_(0, 0) {}

    /**
     * Initialise the list
     * @param frame The current frame
     * @param image The image values
     * @param mask The pixel mask
     */
    PixelList(int frame,
              const af::const_ref<int, af::c_grid<2> > &image,
              const af::const_ref<bool, af::c_grid<2> > &mask) {
      DIALS_ASSERT(image.accessor().all_eq(mask.accessor()));

      frame_ = frame;
      size_ = image.accessor();

      std::size_t num = 0;
      for (std::size_t i = 0; i < mask.size(); ++i) {
        if (mask[i]) num++;
      }

      value_.resize(num);
      index_.resize(num);

      for (std::size_t i = 0, j = 0; i < mask.size(); ++i) {
        if (mask[i]) {
          value_[j] = image[i];
          index_[j] = i;
          j++;
        }
      }
    }

    /**
     * Initialise the list
     * @param frame The current frame
     * @param image The image values
     * @param mask The pixel mask
     */
    PixelList(int frame,
              const af::const_ref<double, af::c_grid<2> > &image,
              const af::const_ref<bool, af::c_grid<2> > &mask) {
      DIALS_ASSERT(image.accessor().all_eq(mask.accessor()));

      frame_ = frame;
      size_ = image.accessor();

      std::size_t num = 0;
      for (std::size_t i = 0; i < mask.size(); ++i) {
        if (mask[i]) num++;
      }

      value_.resize(num);
      index_.resize(num);

      for (std::size_t i = 0, j = 0; i < mask.size(); ++i) {
        if (mask[i]) {
          value_[j] = image[i];
          index_[j] = i;
          j++;
        }
      }
    }

    /**
     * Initialise the list
     * @param frame The current frame
     * @param size The size of the image
     * @param value The list of pixel values
     * @param index The list of pixel indices
     */
    PixelList(int frame,
              int2 size,
              const af::const_ref<double> value,
              const af::const_ref<std::size_t> index)
        : frame_(frame),
          size_(size),
          value_(value.begin(), value.end()),
          index_(index.begin(), index.end()) {
      DIALS_ASSERT(value.size() == index.size());
    }

    /**
     * @return the current frame
     */
    int frame() const {
      return frame_;
    }

    /**
     * @return the image size
     */
    int2 size() const {
      return size_;
    }

    /**
     * @return the pixel values
     */
    af::shared<double> value() const {
      return value_;
    }

    /**
     * @return the pixel indices
     */
    af::shared<std::size_t> index() const {
      return index_;
    }

    /**
     * @return the number of points
     */
    std::size_t num_points() const {
      DIALS_ASSERT(value_.size() == index_.size());
      return value_.size();
    }
  };

  /**
   * A class to label the pixels as spots
   */
  class PixelListLabeller {
  public:
    PixelListLabeller() : size_(0, 0), first_frame_(0), last_frame_(0) {}

    /**
     * Add a pixel list
     * @param pixel_list The pixel list
     */
    void add(const PixelList &pixel_list) {
      // Check the frame number
      if (last_frame_ == first_frame_) {
        first_frame_ = pixel_list.frame();
        size_ = pixel_list.size();
      } else {
        DIALS_ASSERT(pixel_list.frame() == last_frame_);
        DIALS_ASSERT(pixel_list.size().all_eq(size_));
      }

      // Update the last frame number
      int frame = pixel_list.frame();
      last_frame_ = frame + 1;

      // Get pixel list info
      int2 size = pixel_list.size();
      af::const_ref<double> value = pixel_list.value().const_ref();
      af::const_ref<std::size_t> index = pixel_list.index().const_ref();
      DIALS_ASSERT(value.size() == index.size());

      // Add the values and coordinates
      for (std::size_t i = 0; i < value.size(); ++i) {
        std::size_t k = index[i];
        DIALS_ASSERT(i < size[0] * size[1]);
        std::size_t y = k / size[1];
        std::size_t x = k - y * size[1];
        DIALS_ASSERT(x < size[1]);
        DIALS_ASSERT(y < size[0]);
        vec3<int> coord(frame, y, x);
        values_.push_back(value[i]);
        coords_.push_back(coord);
      }
    }

    /**
     * @returns The image size
     */
    int2 size() const {
      return size_;
    }

    /** @returns The number of pixels */
    std::size_t num_pixels() const {
      return values_.size();
    }

    /** @returns The first frame number */
    int first_frame() const {
      return first_frame_;
    }

    /** @returns The last frame number */
    int last_frame() const {
      return last_frame_;
    }

    /** @returns The frame range */
    int2 frame_range() const {
      DIALS_ASSERT(last_frame_ >= first_frame_);
      return int2(first_frame_, last_frame_);
    }

    /** @returns The number of frames */
    std::size_t num_frames() const {
      DIALS_ASSERT(last_frame_ >= first_frame_);
      return last_frame_ - first_frame_;
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
    af::shared<double> values() const {
      return values_;
    }

    /**
     * Label the pixels in 3D
     */
    af::shared<int> labels_3d() const {
      // Adjacency list type
      typedef boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS>
        AdjacencyList;

      if (coords_.size() == 0) {
        return af::shared<int>();
      }

      // Calculate the coordinate indices
      for (std::size_t i = 0; i < coords_.size(); ++i) {
        if (i > 0) {
          DIALS_ASSERT(detail::lessthan(coords_[i - 1], coords_[i]));
        }
      }

      // Create a graph of coordinates
      AdjacencyList graph(coords_.size());
      std::size_t i1 = 0, i2 = 0, i3 = 0;
      for (; i1 < coords_.size() - 1; ++i1) {
        vec3<int> a0 = coords_[i1];
        vec3<int> a1(a0[0], a0[1], a0[2] + 1);
        vec3<int> a2(a0[0], a0[1] + 1, a0[2]);
        vec3<int> a3(a0[0] + 1, a0[1], a0[2]);
        if (coords_[i1 + 1] == a1) {
          boost::add_edge(i1, i1 + 1, graph);
        }
        if (a0[1] < size_[0] - 1) {
          for (; detail::lessthan(coords_[i2], a2) && i2 < coords_.size() - 1; ++i2)
            ;
          if (coords_[i2] == a2) {
            boost::add_edge(i1, i2, graph);
          }
        }
        if (a0[0] < last_frame_ - 1) {
          if (i2 > i3) i3 = i2;
          for (; detail::lessthan(coords_[i3], a3) && i3 < coords_.size() - 1; ++i3)
            ;
          if (coords_[i3] == a3) {
            boost::add_edge(i1, i3, graph);
          }
        }
      }

      // Do the connected components
      af::shared<int> labels(num_vertices(graph), af::init_functor_null<int>());
      DIALS_ASSERT(labels.size() == coords_.size());
      int num = boost::connected_components(graph, &labels[0]);
      DIALS_ASSERT(num <= labels.size());
      return labels;
    }

    /**
     * Label the pixels in 2D
     */
    af::shared<int> labels_2d() const {
      // Adjacency list type
      typedef boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS>
        AdjacencyList;

      if (coords_.size() == 0) {
        return af::shared<int>();
      }

      // Calculate the coordinate indices
      for (std::size_t i = 0; i < coords_.size(); ++i) {
        if (i > 0) {
          DIALS_ASSERT(detail::lessthan(coords_[i - 1], coords_[i]));
        }
      }

      // Create a graph of coordinates
      AdjacencyList graph(coords_.size());
      std::size_t i1 = 0, i2 = 0;
      for (; i1 < coords_.size() - 1; ++i1) {
        vec3<int> a0 = coords_[i1];
        vec3<int> a1(a0[0], a0[1], a0[2] + 1);
        vec3<int> a2(a0[0], a0[1] + 1, a0[2]);
        if (coords_[i1 + 1] == a1) {
          boost::add_edge(i1, i1 + 1, graph);
        }
        if (a0[1] < size_[0] - 1) {
          for (; detail::lessthan(coords_[i2], a2) && i2 < coords_.size() - 1; ++i2)
            ;
          if (coords_[i2] == a2) {
            boost::add_edge(i1, i2, graph);
          }
        }
      }

      // Do the connected components
      af::shared<int> labels(num_vertices(graph));
      DIALS_ASSERT(labels.size() == coords_.size());
      int num = boost::connected_components(graph, &labels[0]);
      DIALS_ASSERT(num <= labels.size());
      return labels;
    }

  private:
    int2 size_;
    int first_frame_;
    int last_frame_;
    af::shared<vec3<int> > coords_;
    af::shared<double> values_;
  };

}}  // namespace dials::model

#endif /* DIALS_MODEL_DATA_PIXEL_LIST_H */
