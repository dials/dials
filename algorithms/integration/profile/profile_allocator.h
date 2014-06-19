/*
 * profile_allocator.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */

#ifndef DIALS_ALGORITHMS_INTEGRATION_PROFILE_PROFILE_ALLOCATOR_H
#define DIALS_ALGORITHMS_INTEGRATION_PROFILE_PROFILE_ALLOCATOR_H

#include <stack>
#include <vector>
#include <scitbx/array_family/tiny_types.h>
#include <scitbx/vec3.h>
#include <cctbx/miller.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  using scitbx::af::int4;
  using scitbx::af::int6;
  using scitbx::vec3;

  /**
   * Class to handle the allocation of profiles with efficient memory usage.
   * This class takes a load of frames and bounding boxes. The algorithm checks
   * which profiles are need for each frame based on a radius for profile
   * fitting and does a single allocation based on the maximum size needed.
   */
  class ProfileAllocator {
  public:

    /**
     * Sort by Z
     */
    struct sort_by_z {
      af::const_ref<int> frame;
      sort_by_z(const af::const_ref<int> &frame_)
        : frame(frame_) {}
      template <typename T>
      bool operator()(T a, T b) {
        return frame[a] < frame[b];
      }
    };

    /**
     * Initialise the allocator.
     */
    ProfileAllocator()
        : max_size_(0),
          max_num_(0) {}

    /**
     * Initialise and allocate memory
     * @param frame The frame numbers
     * @param bbox The bounding boxes
     * @param radius The radius of reflections needed.
     */
    ProfileAllocator(
          const af::const_ref<int> &frame,
          const af::const_ref<int4> &bbox,
          std::size_t radius)
        : max_size_(0),
          max_num_(0),
          max_image_(0),
          offset_(frame.size()) {

      DIALS_ASSERT(frame.size() == bbox.size());

      // Compute maximum size of the 2D shoebox
      for (std::size_t i = 0; i < bbox.size(); ++i) {
        const int4 &b = bbox[i];
        std::size_t size = (b[3] - b[2])*(b[1] - b[0]);
        if (size > max_size_) {
          max_size_ = size;
        }
      }

      // Index of shoebox sorted by min and max z
      std::vector<std::size_t> index(frame.size());
      for (std::size_t i = 0; i < index.size(); ++i) {
        index[i] = i;
      }
      std::sort(index.begin(), index.end(), sort_by_z(frame));

      // Find a position available in the buffer and assign to shoebox. If no
      // position is available add another. When positions become available add
      // them back to the stack of available positions
      std::stack<std::size_t> available;
      std::size_t imin = 0;
      std::size_t imax = 0;
      int length = 2 * radius + 1;
      while (imax < index.size()) {
        int f1 = frame[index[imax]];
        int f0 = f1 - length;
        while (imin < index.size() && frame[index[imin]] <= f0) {
          available.push(offset_[index[imin]]);
          imin++;
        }
        std::size_t count = 0;
        while (imax < index.size() && frame[index[imax]] == f1) {
          std::size_t pos = 0;
          if (available.empty()) {
            pos = max_num_++;
          } else {
            pos = available.top();
            available.pop();
          }
          offset_[index[imax]] = pos;
          imax++;
          count++;
        }
        max_image_ = std::max(max_image_, count);
      }

      // Allocate memory for array
      data_.resize(max_num_ * max_size_, 0);
      lock_.resize(max_num_, -1);
    }

    /**
     * Hold the selected profile.
     */
    void hold(std::size_t index) {
      DIALS_ASSERT(index < offset_.size());
      DIALS_ASSERT(lock_[offset_[index]] == -1);
      lock_[offset_[index]] = index;
    }

    /**
     * Free the selected profile.
     */
    void free(std::size_t index) {
      DIALS_ASSERT(index < offset_.size());
      DIALS_ASSERT(lock_[offset_[index]] == index);
      lock_[offset_[index]] = -1;
    }

    /**
     * Get the profile for the selected index.
     */
    int* data(std::size_t index) {
      DIALS_ASSERT(index < offset_.size());
      DIALS_ASSERT(lock_[offset_[index]] == index);
      return &data_[offset_[index] * max_size_];
    }

    /**
     * Check that the profile is held
     */
    bool held(std::size_t index) {
      DIALS_ASSERT(index < offset_.size());
      return lock_[offset_[index]] == index;
    }

    /**
     * @returns The maximum profile size.
     */
    std::size_t max_size() const {
      return max_size_;
    }

    /**
     * @returns The maximum number of profiles at one time.
     */
    std::size_t max_num() const {
      return max_num_;
    }

    /**
     * @returns The maximum number of profiles in one image.
     */
    std::size_t max_image() const {
      return max_image_;
    }

  private:

    std::size_t max_size_;
    std::size_t max_num_;
    std::size_t max_image_;
    std::vector<int> offset_;
    std::vector<int> lock_;
    std::vector<int> data_;
  };

}}

#endif // DIALS_ALGORITHMS_INTEGRATION_PROFILE_PROFILE_ALLOCATOR_H
