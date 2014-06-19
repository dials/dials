

#ifndef DIALS_ALGORITHMS_INTEGRATION_INTEGRATOR_2D_H
#define DIALS_ALGORITHMS_INTEGRATION_INTEGRATOR_2D_H

#include <stack>
#include <vector>
#include <scitbx/array_family/tiny_types.h>
#include <cctbx/miller.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  using scitbx::af::int4;

  class ProfileAllocator {
  public:

    struct sort_by_z {
      af::const_ref<int> frame;
      sort_by_z(const af::const_ref<int> &frame_)
        : frame(frame_) {}
      template <typename T>
      bool operator()(T a, T b) {
        return frame[a] < frame[b];
      }
    };
    
    ProfileAllocator(
          const af::const_ref<int> &frame,
          const af::const_ref<int4> &bbox,
          double radius)
        : max_size_(0),
          max_num_(0),
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
      int f1 = frame[index[imax]];
      while (true) {
        while (imin < index.size() && frame[index[imin]] <= f1 - 2*radius+1) {
          available.push(offset_[index[imin]]);
          imin++;
        }
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
        }
      }

      // Allocate memory for array
      data_.resize(max_num_ * max_size_, 0);
      lock_.resize(max_num_, -1);
    }

    void hold(std::size_t index) {
      DIALS_ASSERT(index < offset_.size());
      DIALS_ASSERT(lock_[offset_[index]] == -1);
      lock_[offset_[index]] = index;
    }

    void free(std::size_t index) {
      DIALS_ASSERT(index < offset_.size());
      DIALS_ASSERT(lock_[offset_[index]] == index);
      lock_[offset_[index]] = -1;
    }

    int* data(std::size_t index) {
      DIALS_ASSERT(index < offset_.size());
      DIALS_ASSERT(lock_[offset_[index]] == index);
      return &data_[offset_[index]];
    }

  private:

    std::size_t max_size_;
    std::size_t max_num_;
    std::vector<int> offset_;
    std::vector<int> lock_;
    std::vector<int> data_;
  };

  class Integrator2DSpec {
  public:
    
  };
  
  class Integrator2D {
  public:

    Integrator2D(const Integrator2DSpec &spec) {
    }



  };
}}

#endif // DIALS_ALGORITHMS_INTEGRATION_INTEGRATOR_2D_H
