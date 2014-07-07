

#ifndef DIALS_ALGORITHMS_INTEGRATION_INTEGRATOR_2D_H
#define DIALS_ALGORITHMS_INTEGRATION_INTEGRATOR_2D_H

#include <stack>
#include <vector>
#include <scitbx/array_family/tiny_types.h>
#include <scitbx/vec3.h>
#include <cctbx/miller.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/algorithms/integration/profile/profile_allocator.h>
//#include <dials/algorithms/integration/profile/profile_manager.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  using scitbx::af::int4;
  using scitbx::af::int6;
  using scitbx::vec3;

  class Integrator2DSpec {
  public:

    // Data needed for integration
    af::shared< std::size_t > panel;
    af::shared< vec3<double> > xyz;
    af::shared< int6 > bbox;

    // Integration parameters
    std::size_t radius_z;
    std::size_t radius_xy;

    Integrator2DSpec()
      : radius_z(0),
        radius_xy(0) {}

    bool is_valid() const {
      return
        panel.size() > 0 &&
        panel.size() == xyz.size() &&
        panel.size() == bbox.size() &&
        radius_xy > 0;
    }
  };

  class Integrator2D {
  public:

    Integrator2D(const Integrator2DSpec &spec)
      : spec_(spec) {

      // Ensure the spec is valid
      DIALS_ASSERT(spec_.is_valid());

      // Compute the number of partial reflections
      std::size_t size = 0;
      for (std::size_t i = 0; i < spec_.bbox.size(); ++i) {
        int z0 = spec_.bbox[i][4];
        int z1 = spec_.bbox[i][5];
        DIALS_ASSERT(z1 > z0);
        size += z1 - z0;
      }

      // Set the reflection indices, frames and 2d bboxes
      indices_.resize(size);
      af::shared<int> frames(size);
      af::shared<int4> bbox2d(size);
      for (std::size_t i = 0, j = 0; i < spec_.bbox.size(); ++i) {
        int z0 = spec_.bbox[i][4];
        int z1 = spec_.bbox[i][5];
        int4 b(spec_.bbox[i][0], spec_.bbox[i][1],
               spec_.bbox[i][2], spec_.bbox[i][3]);
        for (std::size_t z = z0; z < z1; ++z, ++j) {
          indices_[j] = i;
          frames[j] = z;
          bbox2d[j] = b;
        }
      }

      // Initialise the profile allocator
      profiles_ = ProfileAllocator(
          frames.const_ref(),
          bbox2d.const_ref(),
          spec_.radius_z);
    }

    //void next(const af::const_ref< int, af::c_grid<2> > &image) {

      //// Iterators for indices
      //ProfileManager::iterator first;
      //ProfileManager::iterator last;

      //// Accumulate image pixels into profiles
      //profiles_.accumulate(image);

      //// Loop through all current reflections and compute their backgrounds and
      //// summed intensities. This only requires the lone shoebox to perform.
      //first = profile_.current_begin();
      //last = profile_.current_end();
      //for (ProfileManager::iterator it = first; it != last; ++it) {
        //compute_background(*it);
        //compute_summed_intensity(*it);
      //}

      //// Loop through all reflections that are ready to be profile fitted. This
      //// requires knowledge of reflections on a few adjacent frames
      //first = profile_.ready_begin();
      //last = profile_.ready_end();
      //for (ProfileManager::iterator it = first; it != last; ++it) {
        //compute_profile_fitted_intensity(*it);
      //}

      //// Tell the profile manager to update the frame number
      //profiles_.next();
    //}

    //bool finished() const {
      //return profiles_.finished();
    //}

  private:

    void compute_background(std::size_t index) {

    }

    void compute_summed_intensity(std::size_t index) {

    }

    Integrator2DSpec spec_;
    af::shared<std::size_t> indices_;
    ProfileAllocator profiles_;
  };
}}

#endif // DIALS_ALGORITHMS_INTEGRATION_INTEGRATOR_2D_H
