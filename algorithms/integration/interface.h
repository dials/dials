/*
 * interface.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_INTEGRATION_INTERFACE_H
#define DIALS_ALGORITHMS_INTEGRATION_INTERFACE_H

#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cctbx/miller/index_generator.h>
#include <cctbx/uctbx.h>
#include <cctbx/sgtbx/space_group.h>
#include <dials/array_family/reflection_table.h>
#include <dials/algorithms/reflection_basis/coordinate_system.h>


namespace dials { namespace algorithms {

  class IntegrationManagerData3D {
  public:

    IntegrationManagerData3D(
        af::reflection_table reflections,
        vec2<double> oscillation,
        vec2<int> array_range,
        double block_size) {
      compute_blocks(oscillation, array_range, block_size);
      finished_.assign(blocks_.size(), false);
    }

    af::reflection_table data() {
      DIALS_ASSERT(finished());
      return data_;
    }

    bool finished() const {
      return finished_.all_eq(true);
    }

    std::size_t size() const {
      return finished_.size();
    }

    vec2<int> block(std::size_t index) const {
      DIALS_ASSERT(index < blocks_.size());
      return blocks_[index];
    }

    af::reflection_table operator[](std::size_t index) {
      DIALS_ASSERT(index < finished_.size());
      af::reflection_table result;

      return result;
    }

    void accumulate(std::size_t index, af::reflection_table) {
      DIALS_ASSERT(index < finished_.size());
      DIALS_ASSERT(finished_[index] == false);

      finished_[index] = true;
    }

  private:

    void compute_blocks(
        vec2<double> oscillation,
        vec2<int> array_range,
        double block_size) {
      double phi0 = oscillation[0];
      double dphi = oscillation[1];
      int frame0 = array_range[0];
      int frame1 = array_range[1];
      DIALS_ASSERT(frame1 > frame0);
      int nframes = frame1 - frame0;
      double half_block_size = block_size / 2.0;
      DIALS_ASSERT(half_block_size >= std::abs(dphi));
      DIALS_ASSERT(half_block_size <= std::abs(nframes * dphi));
      double half_block_length_f = half_block_size / dphi;
      int nblocks = (int)std::ceil(nframes / half_block_length_f);
      DIALS_ASSERT(nblocks > 0 && nblocks <= nframes);
      int half_block_length = (int)std::ceil((double)nframes / (double)nblocks);
      af::shared<int> indices;
      indices.push_back(frame0);
      for (int i = 0; i < nblocks; ++i) {
        int frame = frame0 + (i + 1) * half_block_length;
        if (frame > frame1) {
          frame = frame1;
        }
        indices.push_back(frame);
        if (frame == frame1) {
          break;
        }
      }
      DIALS_ASSERT(indices.front() == frame0);
      DIALS_ASSERT(indices.back() == frame1);
      DIALS_ASSERT(indices.size() > 2);
      for (std::size_t i = 0; i < indices.size() - 2; ++i) {
        blocks_.push_back(vec2<int>(indices[i], indices[i+2]));
      }
    }

    af::shared< vec2<int> > blocks_;
    af::shared<bool> finished_;
    af::reflection_table data_;
  };

}}

#endif // DIALS_ALGORITHMS_INTEGRATION_INTERFACE_H
