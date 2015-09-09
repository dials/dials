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

  using scitbx::af::int2;
  using scitbx::vec3;

  /**
   * A class to hold a load of pixels and coordinates
   */
  class PixelList {
  public:

    PixelList() :
       size_(0, 0),
       first_frame_(0) {}

    /**
     * Initialise
     * @param size The size of the 2D image
     * @param first_frame The index of the first frame
     */
    PixelList(int2 size, int first_frame)
      : size_(size),
        first_frame_(first_frame),
        last_frame_(first_frame) {}

    /**
     * Initialise
     * @param size The size of the 2D image
     * @param frame_range The range of the frame indices
     * @param values The image values
     * @param coords The image coords
     */
    PixelList(int2 size, int2 frame_range,
              af::shared<double> values,
              af::shared< vec3<int> > coords)
      : size_(size),
        first_frame_(frame_range[0]),
        last_frame_(frame_range[1]),
        coords_(coords),
        values_(values) {
      DIALS_ASSERT(size[0] > 0 && size[1] > 0);
      DIALS_ASSERT(last_frame_ >= first_frame_);
      DIALS_ASSERT(coords_.size() == values_.size());
      for (std::size_t i = 0; i < coords_.size(); ++i) {
        DIALS_ASSERT(coords[i][0] >= first_frame_ && coords[i][0] < last_frame_);
        DIALS_ASSERT(coords[i][1] >= 0 && coords[i][1] < size_[0]);
        DIALS_ASSERT(coords[i][2] >= 0 && coords[i][2] < size_[1]);
      }
    }

    /** @returns The size of the 2D image */
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
     * Add an image
     * @param image The image pixels
     * @param mask The mask values
     */
    void add_int_image(const af::const_ref< int, af::c_grid<2> > &image,
                       const af::const_ref< bool, af::c_grid<2> > &mask) {
      DIALS_ASSERT(image.accessor().all_eq(mask.accessor()));
      DIALS_ASSERT(image.accessor().all_eq(size_));
      for (std::size_t j = 0; j < size_[0]; ++j) {
        for (std::size_t i = 0; i < size_[1]; ++i) {
          if (mask(j, i)) {
            coords_.push_back(vec3<int>(last_frame_, j, i));
            values_.push_back((double)image(j, i));
          }
        }
      }
      last_frame_++;
    }

    /**
     * Add an image
     * @param image The image pixels
     * @param mask The mask values
     */
    void add_double_image(const af::const_ref< double, af::c_grid<2> > &image,
                          const af::const_ref< bool, af::c_grid<2> > &mask) {
      DIALS_ASSERT(image.accessor().all_eq(mask.accessor()));
      DIALS_ASSERT(image.accessor().all_eq(size_));
      for (std::size_t j = 0; j < size_[0]; ++j) {
        for (std::size_t i = 0; i < size_[1]; ++i) {
          if (mask(j, i)) {
            coords_.push_back(vec3<int>(last_frame_, j, i));
            values_.push_back(image(j, i));
          }
        }
      }
      last_frame_++;
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
    af::shared<double> values() const {
      return values_;
    }

    af::shared<int> labels_3d() const {

      // Adjacency list type
      typedef boost::adjacency_list<
        boost::listS,
        boost::vecS,
        boost::undirectedS> AdjacencyList;

      if (coords_.size() == 0) {
        return af::shared<int>();
      }

      // Calculate the coordinate indices
      std::vector<std::size_t> idx(coords_.size());
      for (std::size_t i = 0; i < coords_.size(); ++i) {
        idx[i] = index(coords_[i]);
        if (i > 0) DIALS_ASSERT(idx[i] > idx[i-1]);
      }

      // Create a graph of coordinates
      AdjacencyList graph(coords_.size());
      std::size_t i1 = 0, i2 = 0, i3 = 0;
      for (; i1 < coords_.size()-1; ++i1) {
        std::size_t idx0 = idx[i1];
        std::size_t idx1 = idx0+1;
        std::size_t idx2 = idx0+size_[1];
        std::size_t idx3 = idx0+size_[0]*size_[1];
        if (idx[i1+1] == idx1 && coords_[i1][2] < size_[1]-1) {
          boost::add_edge(i1, i1+1, graph);
        }
        if (coords_[i1][1] < size_[0]-1) {
          for (; idx[i2] < idx2 && i2 < coords_.size(); ++i2);
          if (idx[i2] == idx2) {
            boost::add_edge(i1, i2, graph);
          }
        }
        if (coords_[i1][0] < last_frame_-1) {
          if (i2 > i3) i3 = i2;
          for (; idx[i3] < idx3 && i3 < coords_.size(); ++i3);
          if (idx[i3] == idx3) {
            boost::add_edge(i1, i3, graph);
          }
        }
      }

      // Do the connected components
      af::shared<int> labels(num_vertices(graph), af::init_functor_null<int>());
      int num = boost::connected_components(graph, &labels[0]);
      DIALS_ASSERT(num <= labels.size());
      return labels;
    }

    af::shared<int> labels_2d() const {

      // Adjacency list type
      typedef boost::adjacency_list<
        boost::listS,
        boost::vecS,
        boost::undirectedS> AdjacencyList;

      if (coords_.size() == 0) {
        return af::shared<int>();
      }

      // Calculate the coordinate indices
      std::vector<std::size_t> idx(coords_.size());
      for (std::size_t i = 0; i < coords_.size(); ++i) {
        idx[i] = index(coords_[i]);
        if (i > 0) DIALS_ASSERT(idx[i] > idx[i-1]);
      }

      // Create a graph of coordinates
      AdjacencyList graph(coords_.size());
      std::size_t i1 = 0, i2 = 0;
      for (; i1 < coords_.size()-1; ++i1) {
        std::size_t idx0 = idx[i1];
        std::size_t idx1 = idx0+1;
        std::size_t idx2 = idx0+size_[1];
        if (idx[i1+1] == idx1 && coords_[i1][2] < size_[1]-1) {
          boost::add_edge(i1, i1+1, graph);
        }
        if (coords_[i1][1] < size_[0]-1) {
          for (; idx[i2] < idx2 && i2 < coords_.size(); ++i2);
          if (idx[i2] == idx2) {
            boost::add_edge(i1, i2, graph);
          }
        }
      }

      // Do the connected components
      af::shared<int> labels(num_vertices(graph));
      int num = boost::connected_components(graph, &labels[0]);
      DIALS_ASSERT(num <= labels.size());
      return labels;
    }

  private:

    std::size_t index(vec3<int> a) const {
      DIALS_ASSERT(a[0] >= first_frame_ && a[0] < last_frame_);
      DIALS_ASSERT(a[1] >= 0 && a[1] < size_[0]);
      DIALS_ASSERT(a[2] >= 0 && a[2] < size_[1]);
      return ((std::size_t)(a[0] - first_frame_) * (std::size_t)size_[0] +
              (std::size_t)a[1]) * (std::size_t)size_[1] +
              (std::size_t)a[2];
    }

    int2 size_;
    int first_frame_;
    int last_frame_;
    af::shared< vec3<int> > coords_;
    af::shared<double> values_;
  };

}} // namespace dials::model

#endif /* DIALS_MODEL_DATA_PIXEL_LIST_H */
