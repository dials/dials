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
#include <numeric>
#include <list>
#include <vector>
#include <dials/model/data/image.h>
#include <dials/model/data/shoebox.h>
#include <dials/array_family/reflection_table.h>
#include <dials/array_family/boost_python/flex_table_suite.h>


namespace dials { namespace algorithms {

  using model::Image;
  using model::Shoebox;

  /**
   * A class to compute integration jobs.
   */
  class IntegrationJobCalculator {
  public:

    /**
     * Compute the integration jobs
     * @param array_range The range of frames to process
     * @param block_size The number of frames in a job
     */
    IntegrationJobCalculator(
        vec2<int> array_range,
        double block_size) {
      int frame0 = array_range[0];
      int frame1 = array_range[1];
      DIALS_ASSERT(frame1 > frame0);
      int nframes = frame1 - frame0;
      DIALS_ASSERT(nframes > 0);
      if (block_size > nframes) {
        block_size = nframes;
      }
      DIALS_ASSERT(block_size > 0);
      if (block_size == 1) {
        for (int f = frame0; f < frame1; ++f) {
          jobs_.push_back(tiny<int,2>(f, f+1));
        }
      } else {
        int nblocks = (int)std::ceil(2.0 * nframes / block_size);
        DIALS_ASSERT(nblocks > 0 && nblocks <= nframes);
        int half_block_size = (int)std::ceil((double)nframes / (double)nblocks);
        af::shared<int> indices;
        indices.push_back(frame0);
        for (int i = 0; i < nblocks; ++i) {
          int frame = frame0 + (i + 1) * half_block_size;
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
          int i1 = indices[i];
          int i2 = indices[i+2];
          DIALS_ASSERT(i2 > i1);
          jobs_.push_back(tiny<int,2>(i1, i2));
        }
        DIALS_ASSERT(jobs_.size() > 0);
      }
    }

    /**
     * @returns The list of jobs.
     */
    af::shared< tiny<int,2> > jobs() const {
      return af::shared< tiny<int,2> >(&jobs_[0], &jobs_[0] + jobs_.size());
    }

  private:
    std::vector< tiny<int,2> > jobs_;
  };


  /**
   * A class to managing spliting and mergin data
   */
  class IntegrationManagerExecutor {
  public:

    IntegrationManagerExecutor(
        const IntegrationJobCalculator &jobcalculator,
        af::reflection_table reflections)
          : jobs_(jobcalculator.jobs()),
            data_(reflections) {
      DIALS_ASSERT(jobs_.size() > 0);

      // Set all the finished flags to false
      finished_.assign(jobs_.size(), false);

      // Generate indices of reflections to be integrated, used as reference
      // spots or passed just as data for each data block. If the reflection is
      // not to be integrated, it is added to each block which it overlaps. If
      // the reflection is a reference spot, it is added to each block in which
      // it is fully recorded. If the spot is to be integrated, it is added to
      // the block in which it is closest to the centre. If the reflection is
      // larger than block_size / 2, then it is not fully recorded in any block
      // and is unprocessed.
      typedef std::pair<std::size_t, bool> job_type;
      typedef std::vector<job_type> job_list_type;

      // Get some reflection data
      af::const_ref<int6> bbox = data_["bbox"];
      af::ref<std::size_t> flags = data_["flags"];

      // Get which reflections to process in which job and task
      std::vector<job_list_type> indices(jobs_.size());
      for (std::size_t index = 0; index < bbox.size(); ++index) {
        int z0 = bbox[index][4];
        int z1 = bbox[index][5];
        std::size_t &f = flags[index];
        if (!(f & af::DontIntegrate)) {
          std::size_t jmin = 0;
          double dmin = 0;
          bool first = true;
          for (std::size_t j = 0; j < jobs_.size(); ++j) {
            int jz0 = jobs_[j][0];
            int jz1 = jobs_[j][1];
            if (z0 >= jz0 && z1 <= jz1) {
              if (f & af::ReferenceSpot) {
                indices[j].push_back(job_type(index, false));
              }
              double zc = (z1 + z0) / 2.0;
              double jc = (jz1 + jz0) / 2.0;
              double d = std::abs(zc - jc);
              if (first || d < dmin) {
                jmin = j;
                dmin = d;
                first = false;
              }
            }
          }
          int jz0 = jobs_[jmin][0];
          int jz1 = jobs_[jmin][1];
          if (first == false && z0 >= jz0 && z1 <= jz1) {
            if (f & af::ReferenceSpot) {
              DIALS_ASSERT(indices[jmin].size() > 0);
              DIALS_ASSERT(indices[jmin].back().first == index);
              indices[jmin].back().second = true;
            } else {
              indices[jmin].push_back(job_type(index, true));
            }
          } else {
            f |= af::DontIntegrate;
            ignored_.push_back(index);
          }
        } else {
          ignored_.push_back(index);
        }
      }

      // Compute number of reflections in each task
      std::vector<std::size_t> num(jobs_.size(), 0);
      for (std::size_t i = 0; i < num.size(); ++i) {
        num[i] = indices[i].size();
      }

      // Compute offsests
      offset_.push_back(0);
      std::partial_sum(num.begin(), num.end(), std::back_inserter(offset_));

      // Compute indices
      indices_.resize(offset_.back());
      mask_.resize(offset_.back());
      std::size_t k = 0;
      for (std::size_t i = 0; i < indices.size(); ++i) {
        const job_list_type& ind = indices[i];
        for (std::size_t j = 0; j < ind.size(); ++j, ++k) {
          DIALS_ASSERT(k < indices_.size());
          indices_[k] = ind[j].first;
          mask_[k] = ind[j].second;
          DIALS_ASSERT(indices_[k] < bbox.size());
        }
      }
      DIALS_ASSERT(k == indices_.size());
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
    vec2<int> job(std::size_t index) const {
      DIALS_ASSERT(index < jobs_.size());
      return jobs_[index];
    }

    /**
     * @returns The list of reflections to not process.
     */
    af::shared<std::size_t> ignored() const {
      return af::shared<std::size_t>(&ignored_[0], &ignored_[0]+ignored_.size());
    }

    /**
     * @returns The reflections for a particular block.
     */
    af::reflection_table split(std::size_t index) {

      using namespace af::boost_python::flex_table_suite;

      // Check the input
      DIALS_ASSERT(index < finished_.size());
      af::const_ref<std::size_t> ind = indices(index);
      af::const_ref<bool> msk = mask(index);
      DIALS_ASSERT(ind.size() == msk.size());

      // Extract the reflection table
      af::reflection_table result = select_rows_index(data_, ind);

      // Extract the flags and set those reflections that are not to be
      // processed.
      af::ref<std::size_t> flags = result["flags"];
      for (std::size_t i = 0; i < flags.size(); ++i) {
        if (msk[i] == false) {
          flags[i] |= af::DontIntegrate;
        }
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
      af::const_ref<std::size_t> ind = indices(index);
      af::const_ref<bool> msk = mask(index);
      DIALS_ASSERT(ind.size() == msk.size());
      DIALS_ASSERT(ind.size() == result.size());

      // Set the result
      set_selected_rows_index_mask(data_, ind, msk, result);

      // Set finished flag
      finished_[index] = true;
    }

  private:

    /**
     * Get the indices for each job
     */
    af::const_ref<std::size_t> indices(std::size_t index) const {
      DIALS_ASSERT(index < offset_.size()-1);
      std::size_t i0 = offset_[index];
      std::size_t i1 = offset_[index+1];
      DIALS_ASSERT(i1 >= i0);
      std::size_t off = i0;
      std::size_t num = i1 - i0;
      DIALS_ASSERT(off + num <= indices_.size());
      return af::const_ref<std::size_t> (&indices_[off], num);
    }

    /**
     * Get the mask for each job
     */
    af::const_ref<bool> mask(std::size_t index) const {
      DIALS_ASSERT(index < offset_.size()-1);
      std::size_t i0 = offset_[index];
      std::size_t i1 = offset_[index+1];
      DIALS_ASSERT(i1 >= i0);
      std::size_t off = i0;
      std::size_t num = i1 - i0;
      DIALS_ASSERT(off + num <= mask_.size());
      return af::const_ref<bool> (&mask_[off], num);
    }

    af::shared< tiny<int,2> > jobs_;
    af::shared<bool> finished_;
    af::reflection_table data_;
    af::shared<std::size_t> offset_;
    af::shared<std::size_t> indices_;
    af::shared<bool> mask_;
   af::shared<std::size_t> ignored_;
  };

}}

#endif // DIALS_ALGORITHMS_INTEGRATION_INTERFACE_H
