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
   * A class to manage jobs for multiple sweeps
   */
  class JobList {
  public:

    /**
     * A single job class.
     */
    class Job {
    public:

      Job(tiny<int,2> expr, tiny<int,2> frames)
        : expr_(expr),
          frames_(frames) {
        DIALS_ASSERT(expr[0] >= 0 && expr[1] > expr[0]);
        DIALS_ASSERT(frames[1] > frames[0]);
      }

      /** @returns The experiments which this job covers */
      tiny<int,2> expr() const {
        return expr_;
      }

      /** @returns The frames which this job covers */
      tiny<int,2> frames() const {
        return frames_;
      }

      /** @returns The number of experiments which this job covers */
      std::size_t nexpr() const {
        return expr_[1] - expr_[0];
      }

      /** @returns The number of frames which this job covers */
      std::size_t nframes() const {
        return frames_[1] - frames_[0];
      }

    private:
      tiny<int,2> expr_;
      tiny<int,2> frames_;
    };

    JobList() {}

    /**
     * Add a new group of jobs covering a range of experiments
     * @param expr The range of experiments
     * @param range The range of frames
     * @param block_size The job block size
     */
    void add(tiny<int,2> expr, tiny<int,2> range, int block_size) {
      DIALS_ASSERT(expr[0] >= 0);
      DIALS_ASSERT(expr[1] > expr[0]);
      DIALS_ASSERT(expr_.size() == 0 || expr[0] >= expr_.back()[1]);
      add_jobs(expr, range, block_size);
      expr_.push_back(expr);
    }

    /**
     * @returns The requested job
     */
    const Job& operator[](std::size_t index) const {
      DIALS_ASSERT(index < jobs_.size());
      return jobs_[index];
    }

    /**
     * @returns The number of jobs
     */
    std::size_t size() const {
      return jobs_.size();
    }

  private:

    void add_jobs(tiny<int,2> expr, tiny<int,2> range, int block_size) {
      int frame0 = range[0];
      int frame1 = range[1];
      DIALS_ASSERT(frame1 > frame0);
      int nframes = frame1 - frame0;
      DIALS_ASSERT(nframes > 0);
      if (block_size > nframes) {
        block_size = nframes;
      }
      DIALS_ASSERT(block_size > 0);
      if (block_size == 1) {
        for (int f = frame0; f < frame1; ++f) {
          jobs_.push_back(Job(expr, tiny<int,2>(f, f+1)));
        }
      } else {
        int nblocks = (int)std::ceil(2.0 * nframes / (double)block_size);
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
          jobs_.push_back(Job(expr, tiny<int,2>(i1, i2)));
        }
        DIALS_ASSERT(jobs_.size() > 0);
      }
    }

    std::vector<Job> jobs_;
    std::vector< tiny<int,2> > expr_;
  };


  /**
   * A class to managing reflection lookup indices
   */
  class ReflectionLookup {
  public:

    ReflectionLookup(
        const af::const_ref<std::size_t> &id,
        const af::const_ref<std::size_t> &flags,
        const af::const_ref<int6> &bbox,
        const af::const_ref< tiny<int,2> > &jobs)
          : jobs_(jobs.begin(), jobs.end()) {
      DIALS_ASSERT(jobs_.size() > 0);

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

      // Check blocks are valid (can be overlapping)
      DIALS_ASSERT(jobs_[0][1] > jobs_[0][0]);
      for (std::size_t i = 1; i < jobs_.size(); ++i) {
        DIALS_ASSERT(jobs_[i][1] > jobs_[i][0]);
        DIALS_ASSERT(jobs_[i][0] > jobs_[i-1][0]);
        DIALS_ASSERT(jobs_[i][1] > jobs_[i-1][1]);
        DIALS_ASSERT(jobs_[i][0] <= jobs_[i-1][1]);
      }

      // Check all the reflections are in range
      int frame0 = jobs_.front()[0];
      int frame1 = jobs_.back()[1];
      DIALS_ASSERT(frame1 > frame0);
      for (std::size_t i = 0; i < bbox.size(); ++i) {
        DIALS_ASSERT(bbox[i][1] > bbox[i][0]);
        DIALS_ASSERT(bbox[i][3] > bbox[i][2]);
        DIALS_ASSERT(bbox[i][5] > bbox[i][4]);
        DIALS_ASSERT(bbox[i][4] >= frame0);
        DIALS_ASSERT(bbox[i][5] <= frame1);
      }

      // Create the lookups
      std::vector<std::size_t> lookup0(frame1 - frame0);
      std::vector<std::size_t> lookup1(frame1 - frame0);
      int frame = frame0;
      for (std::size_t i = 0; i < jobs_.size(); ++i) {
        tiny<int,2> b = jobs_[i];
        DIALS_ASSERT(frame >= b[0]);
        for (; frame < b[1]; ++frame) {
          lookup0[frame-frame0] = i;
        }
      }
      DIALS_ASSERT(frame == frame1);
      for (std::size_t i = 0; i < jobs_.size(); ++i) {
        std::size_t j = jobs_.size() - i - 1;
        tiny<int,2> b = jobs_[j];
        DIALS_ASSERT(frame <= b[1]);
        for (; frame > b[0]; --frame) {
          lookup1[frame-frame0-1] = j;
        }
      }
      DIALS_ASSERT(frame == frame0);

      // Check the lookups
      for (std::size_t i = 1; i < lookup0.size(); ++i) {
        DIALS_ASSERT(lookup0[i] >= lookup0[i-1]);
        DIALS_ASSERT(lookup1[i] >= lookup1[i-1]);
      }

      // Get which reflections to process in which job and task
      std::vector<job_list_type> indices(jobs_.size());
      for (std::size_t index = 0; index < bbox.size(); ++index) {
        int z0 = bbox[index][4];
        int z1 = bbox[index][5];
        const std::size_t &f = flags[index];
        if (!(f & af::DontIntegrate)) {
          std::size_t j0 = lookup0[z0-frame0];
          std::size_t j1 = lookup1[z1-frame0-1];
          DIALS_ASSERT(j0 < jobs_.size());
          DIALS_ASSERT(j1 < jobs_.size());
          DIALS_ASSERT(j1 >= j0);
          DIALS_ASSERT(z0 >= jobs_[j0][0]);
          DIALS_ASSERT(z1 <= jobs_[j1][1]);
          std::size_t jmin = 0;
          double dmin = 0;
          bool inside = false;
          for (std::size_t j = j0; j <= j1; ++j) {
            int jz0 = jobs_[j][0];
            int jz1 = jobs_[j][1];
            if (z0 >= jz0 && z1 <= jz1) {
              if (f & af::ReferenceSpot) {
                indices[j].push_back(job_type(index, false));
              }
              double zc = (z1 + z0) / 2.0;
              double jc = (jz1 + jz0) / 2.0;
              double d = std::abs(zc - jc);
              if (!inside || d < dmin) {
                jmin = j;
                dmin = d;
                inside = true;
              }
            }
          }
          int jz0 = jobs_[jmin][0];
          int jz1 = jobs_[jmin][1];
          DIALS_ASSERT(inside == true);
          DIALS_ASSERT(z0 >= jz0 && z1 <= jz1);
          if (f & af::ReferenceSpot) {
            DIALS_ASSERT(indices[jmin].size() > 0);
            DIALS_ASSERT(indices[jmin].back().first == index);
            indices[jmin].back().second = true;
          } else {
            indices[jmin].push_back(job_type(index, true));
          }
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
     * @returns The number of tasks
     */
    std::size_t size() const {
      return jobs_.size();
    }

    /**
     * @returns The block indices
     */
    vec2<int> job(std::size_t index) const {
      DIALS_ASSERT(index < jobs_.size());
      return jobs_[index];
    }

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
    af::shared<std::size_t> offset_;
    af::shared<std::size_t> indices_;
    af::shared<bool> mask_;
  };


  /**
   * A class to managing spliting and mergin data
   */
  class ReflectionManager {
  public:

    /**
     * Create the reflection manager
     * @param jobs The job calculator
     * @param groups The group that each experiment is in
     * @param data The reflection data
     */
    ReflectionManager(
        const JobList &jobs,
        af::reflection_table data)
          : lookup_(init(jobs, data)),
            data_(data),
            finished_(lookup_.size(), false) {
      DIALS_ASSERT(finished_.size() > 0);
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
     * @returns The job
     */
    tiny<int,2> job(std::size_t index) const {
      return lookup_.job(index);
    }

    /**
     * @returns The reflections for a particular block.
     */
    af::reflection_table split(std::size_t index) {

      using namespace af::boost_python::flex_table_suite;

      // Check the input
      DIALS_ASSERT(index < finished_.size());
      af::const_ref<std::size_t> ind = lookup_.indices(index);
      af::const_ref<bool> msk = lookup_.mask(index);
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
      af::const_ref<std::size_t> ind = lookup_.indices(index);
      af::const_ref<bool> msk = lookup_.mask(index);
      DIALS_ASSERT(ind.size() == msk.size());
      DIALS_ASSERT(ind.size() == result.size());

      // Set the result
      set_selected_rows_index_mask(data_, ind, msk, result);

      // Set finished flag
      finished_[index] = true;
    }

  private:

    /**
     * Initialise the indexer
     */
    ReflectionLookup init(
        const JobList &jobs,
        af::reflection_table data) {
      DIALS_ASSERT(data.is_consistent());
      DIALS_ASSERT(data.size() > 0);
      DIALS_ASSERT(data.contains("id"));
      DIALS_ASSERT(data.contains("flags"));
      DIALS_ASSERT(data.contains("bbox"));
      DIALS_ASSERT(jobs.size() > 0);
      af::shared< tiny<int,2> > job_items(jobs.size());
      for (std::size_t i = 0; i < jobs.size(); ++i) {
        job_items[i] = jobs[i].frames();
      }
      return ReflectionLookup(
          data["id"],
          data["flags"],
          data["bbox"],
          job_items.const_ref());
    }

    ReflectionLookup lookup_;
    af::reflection_table data_;
    af::shared<bool> finished_;
  };

}}

#endif // DIALS_ALGORITHMS_INTEGRATION_INTERFACE_H
