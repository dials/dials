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
#include <dials/array_family/boost_python/flex_table_suite.h>
#include <dials/algorithms/reflection_basis/coordinate_system.h>


namespace dials { namespace algorithms {

  /**
   * A class to managing spliting and mergin data
   */
  class IntegrationManagerData3D {
  public:

    IntegrationManagerData3D(
        af::reflection_table reflections,
        vec2<double> oscillation,
        vec2<int> array_range,
        double block_size)
          : data_(reflections) {

      // Compute the blocks
      compute_blocks(oscillation, array_range, block_size);
      finished_.assign(blocks_.size(), false);

      // Get the bounding boxes and flags
      af::const_ref<int6> bbox = reflections["bbox"];
      af::ref<std::size_t> flags = reflections["flags"];

      // Generate indices of reflections to be integrated, used as reference
      // spots or passed just as data for each data block. If the reflection is
      // not to be integrated, it is added to each block which it overlaps. If
      // the reflection is a reference spot, it is added to each block in which
      // it is fully recorded. If the spot is to be integrated, it is added to
      // the block in which it is closest to the centre. If the reflection is
      // larger than block_size / 2, then it is not fully recorded in any block
      // and is unprocessed.
      to_process_.resize(blocks_.size());
      to_include_.resize(blocks_.size());
      for (std::size_t i = 0; i < bbox.size(); ++i) {
        int z0 = bbox[i][4];
        int z1 = bbox[i][5];
        std::size_t &f = flags[i];
        if (!(f & af::DontIntegrate)) {
          std::size_t jmin = 0;
          double dmin = 0;
          for (std::size_t j = 0; j < blocks_.size(); ++j) {
            int bz0 = blocks_[j][0];
            int bz1 = blocks_[j][1];
            if (f & af::ReferenceSpot) {
              if (z0 >= bz0 && z1 <= bz1) {
                to_include_[j].push_back(i);
              }
            }
            double zc = (z1 + z0) / 2.0;
            double bc = (bz1 + bz0) / 2.0;
            double d = (zc - bc)*(zc - bc);
            if (j == 0 || d < dmin) {
              jmin = j;
              dmin = d;
            }
          }
          int bz0 = blocks_[jmin][0];
          int bz1 = blocks_[jmin][1];
          if (z0 >= bz0 && z1 <= bz1) {
            to_process_[jmin].push_back(i);
          } else {
            to_not_process_.push_back(i);
            f |= af::DontIntegrate;
          }
        }
        if (f & af::DontIntegrate) {
          for (std::size_t j = 0; j < blocks_.size(); ++j) {
            int bz0 = blocks_[j][0];
            int bz1 = blocks_[j][1];
            if (!(z1 <= bz0 || z0 >= bz1)) {
              to_include_[j].push_back(i);
            }
          }
        }
      }
    }

    /**
     * @returns The result data
     */
    af::reflection_table data() {
      DIALS_ASSERT(finished());
      return data_;
    }

    /**
     * @returns Is the process finished
     */
    bool finished() const {
      return finished_.all_eq(true);
    }

    /**
     * @returns The number of tasks
     */
    std::size_t size() const {
      return finished_.size();
    }

    /**
     * @returns The block indices
     */
    vec2<int> block(std::size_t index) const {
      DIALS_ASSERT(index < blocks_.size());
      return blocks_[index];
    }

    /**
     * @returns The list of reflections to not process.
     */
    af::shared<std::size_t> to_not_process() const {
      return af::shared<std::size_t>(
          &to_not_process_[0],
          &to_not_process_[0] + to_not_process_.size());
    }

    /**
     * @returns The list of reflections to include.
     */
    af::shared<std::size_t> to_include(std::size_t index) const {
      DIALS_ASSERT(index < blocks_.size());
      return af::shared<std::size_t>(
          &to_include_[index][0],
          &to_include_[index][0] + to_include_[index].size());
    }

    /**
     * @returns The list of reflections to process.
     */
    af::shared<std::size_t> to_process(std::size_t index) const {
      DIALS_ASSERT(index < blocks_.size());
      return af::shared<std::size_t>(
          &to_process_[index][0],
          &to_process_[index][0] + to_process_[index].size());
    }

    /**
     * @returns The reflections for a particular block.
     */
    af::reflection_table split(std::size_t index) {

      using namespace af::boost_python::flex_table_suite;

      // Check the input
      DIALS_ASSERT(index < finished_.size());

      // Get the indices of reflections to select
      std::vector<std::size_t> &process = to_process_[index];
      std::vector<std::size_t> &include = to_include_[index];
      af::shared<std::size_t> indices(process.size() + include.size());
      std::copy(process.begin(), process.end(), indices.begin());
      std::copy(include.begin(), include.end(), indices.begin() + process.size());

      // Extract the reflection table
      af::reflection_table result = select_rows_index(
          data_, indices.const_ref());

      // Extract the flags and set those reflections that are not to be
      // processed.
      af::ref<std::size_t> bk_id = result["bk_id"];
      af::ref<std::size_t> flags = result["flags"];
      for (std::size_t i = 0; i < bk_id.size(); ++i) {
        bk_id[i] = i;
      }
      for (std::size_t i = process.size(); i < flags.size(); ++i) {
        flags[i] |= af::DontIntegrate;
      }

      // Return the reflections
      return result;
    }

    /**
     * Accumulate the results.
     */
    void accumulate(std::size_t index, af::reflection_table result) {

      using namespace af::boost_python::flex_table_suite;

      // Check the input
      DIALS_ASSERT(index < finished_.size());
      DIALS_ASSERT(finished_[index] == false);

      // Get the indices of reflections to select
      std::vector<std::size_t> &process = to_process_[index];
      std::vector<std::size_t> &include = to_include_[index];
      af::const_ref<std::size_t> indices(&process[0], process.size());

      // Resize the input to only select those which should have been processed.
      // Check that the book-keeping indices are as expected
      DIALS_ASSERT(process.size() + include.size() == result.size());
      result.resize(process.size());
      af::const_ref<std::size_t> bk_id = result["bk_id"];
      for (std::size_t i = 0; i < bk_id.size(); ++i) {
        DIALS_ASSERT(bk_id[i] == i);
      }

      // Set the result
      set_selected_rows_index(data_, indices, result);

      // Set finished flag
      finished_[index] = true;
    }

  private:

    /**
     * Compute the blocks
     */
    void compute_blocks(
        vec2<double> oscillation,
        vec2<int> array_range,
        double block_size) {
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
      DIALS_ASSERT(blocks_.size() > 0);
    }

    af::shared< vec2<int> > blocks_;
    af::shared<bool> finished_;
    af::reflection_table data_;
    std::vector< std::vector<std::size_t> > to_process_;
    std::vector< std::vector<std::size_t> > to_include_;
    std::vector<std::size_t> to_not_process_;
  };

}}

#endif // DIALS_ALGORITHMS_INTEGRATION_INTERFACE_H
