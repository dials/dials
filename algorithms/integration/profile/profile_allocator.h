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

    typedef const std::size_t* iterator;

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
          max_num_(0),
          max_image_(0),
          radius_(0) {}

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
          radius_(radius),
          index_(frame.size()),
          offset_(frame.size()) {

      DIALS_ASSERT(frame.size() == bbox.size());

      // Compute maximum size of the 2D shoebox
      std::size_t num_frames = 0;
      for (std::size_t i = 0; i < bbox.size(); ++i) {
        const int4 &b = bbox[i];
        std::size_t size = (b[3] - b[2])*(b[1] - b[0]);
        if (size > max_size_) {
          max_size_ = size;
        }
        DIALS_ASSERT(frame[i] >= 0);
        if (frame[i] > num_frames) {
          num_frames = frame[i];
        }
      }
      num_frames++;

      // Index of shoebox sorted by min and max z
      for (std::size_t i = 0; i < index_.size(); ++i) {
        index_[i] = i;
      }
      std::sort(index_.begin(), index_.end(), sort_by_z(frame));

      // Setup the frame offset index
      frame_offset_.resize(num_frames+1);
      frame_offset_[0] = 0;
      int current_frame = frame[index_[0]];
      DIALS_ASSERT(current_frame == 0);
      for (std::size_t i = 1; i < index_.size(); ++i) {
        if (frame[index_[i]] == current_frame) {
          continue;
        } else if (frame[index_[i]] == current_frame+1) {
          current_frame = frame[index_[i]];
          frame_offset_[current_frame] = i;
        } else {
          DIALS_ERROR("No predicted reflections on frame");
        }
      }
      DIALS_ASSERT(current_frame+1 == frame_offset_.size() - 1);
      frame_offset_[current_frame+1] = index_.size();

      // Find a position available in the buffer and assign to shoebox. If no
      // position is available add another. When positions become available add
      // them back to the stack of available positions
      std::stack<std::size_t> available;
      std::size_t imin = 0;
      std::size_t imax = 0;
      int length = 2 * radius + 1;
      while (imax < index_.size()) {
        int f1 = frame[index_[imax]];
        int f0 = f1 - length;
        while (imin < index_.size() && frame[index_[imin]] <= f0) {
          available.push(offset_[index_[imin]]);
          imin++;
        }
        std::size_t count = 0;
        while (imax < index_.size() && frame[index_[imax]] == f1) {
          std::size_t pos = 0;
          if (available.empty()) {
            pos = max_num_++;
          } else {
            pos = available.top();
            available.pop();
          }
          offset_[index_[imax]] = pos;
          imax++;
          count++;
        }
        max_image_ = std::max(max_image_, count);
      }

      // Allocate memory for array
      data_.resize(max_num_ * max_size_, 0);
      mask_.resize(max_num_ * max_size_, 0);
      bgrd_.resize(max_num_ * max_size_, 0);
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
      DIALS_ASSERT(held(index));
      lock_[offset_[index]] = -1;
    }

    /**
     * Loop through all the reflections that have finished and free them. Then
     * loop through all those that are starting and hold them.
     */
    void lock(std::size_t frame) {
      iterator beg, end;
      beg = begin_free(frame);
      end = end_free(frame);
      for (iterator it = beg; it != end; ++it) {
        free(*it);
      }
      beg = begin_hold(frame);
      end = end_hold(frame);
      for (iterator it = beg; it != end; ++it) {
        hold(*it);
      }
    }

    /**
     * Get the profile for the selected index.
     */
    int* data(std::size_t index) {
      DIALS_ASSERT(held(index));
      return &data_[offset_[index] * max_size_];
    }

    /**
     * Get the profile mask for the selected index.
     */
    int* mask(std::size_t index) {
      DIALS_ASSERT(held(index));
      return &mask_[offset_[index] * max_size_];
    }

    /**
     * Get the profile background for the selected index.
     */
    double* background(std::size_t index) {
      DIALS_ASSERT(held(index));
      return &bgrd_[offset_[index] * max_size_];
    }

    /**
     * Get the first index of partial reflections for this frame.
     */
    iterator begin(std::size_t frame) const {
      DIALS_ASSERT(frame < frame_offset_.size());
      return &index_[frame_offset_[frame]];
    }

    /**
     * Get the last index of partial reflections for this frame.
     */
    iterator end(std::size_t frame) const {
      return begin(frame+1);
    }

    /**
     * Get the first index of partial reflections for this frame.
     */
    iterator begin_active(std::size_t frame) const {
      return begin(frame > 2*radius_ ? frame - 2*radius_ : 0);
    }

    /**
     * Get the last index of partial reflections for this frame.
     */
    iterator end_active(std::size_t frame) const {
      return end(frame);
    }

    /**
     * Get the first index of partial reflections for this frame.
     */
    iterator begin_free(std::size_t frame) const {
      return begin_active(frame > 0 ? frame - 1 : 0);
    }

    /**
     * Get the last index of partial reflections for this frame.
     */
    iterator end_free(std::size_t frame) const {
      return begin_free(frame+1);
    }

    /**
     * Get the first index of partial reflections for this frame.
     */
    iterator begin_hold(std::size_t frame) const {
      return begin(frame);
    }

    /**
     * Get the last index of partial reflections for this frame.
     */
    iterator end_hold(std::size_t frame) const {
      return end(frame);
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
    std::size_t radius_;
    std::vector<std::size_t> index_;
    std::vector<std::size_t> frame_offset_;
    std::vector<int> offset_;
    std::vector<int> lock_;
    std::vector<int> data_;
    std::vector<int> mask_;
    std::vector<double> bgrd_;
  };

}}

#endif // DIALS_ALGORITHMS_INTEGRATION_PROFILE_PROFILE_ALLOCATOR_H
