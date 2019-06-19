/*
 * integrator.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_INTEGRATION_INTEGRATOR_H
#define DIALS_ALGORITHMS_INTEGRATION_INTEGRATOR_H

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
   * A class to manage groups of jobs
   */
  class GroupList {
  public:
    /**
     * A group class
     */
    class Group {
    public:
      Group(tiny<int, 2> index, tiny<int, 2> expr, tiny<int, 2> frames)
          : index_(index), expr_(expr), frames_(frames) {
        DIALS_ASSERT(index[0] >= 0 && index[1] > index[0]);
        DIALS_ASSERT(expr[0] >= 0 && expr[1] > expr[0]);
        DIALS_ASSERT(frames[1] > frames[0]);
      }

      /** @returns The jobs indices */
      tiny<int, 2> index() const {
        return index_;
      }

      /** @returns The experiments which this job covers */
      tiny<int, 2> expr() const {
        return expr_;
      }

      /** @returns The frames which this job covers */
      tiny<int, 2> frames() const {
        return frames_;
      }

      /** @returns The number of jobs which this group covers */
      std::size_t nindex() const {
        return index_[1] - index_[0];
      }

      /** @returns The number of experiments which this group covers */
      std::size_t nexpr() const {
        return expr_[1] - expr_[0];
      }

      /** @returns The number of frames which this group covers */
      std::size_t nframes() const {
        return frames_[1] - frames_[0];
      }

    private:
      tiny<int, 2> index_;
      tiny<int, 2> expr_;
      tiny<int, 2> frames_;
    };

    void add(tiny<int, 2> index, tiny<int, 2> expr, tiny<int, 2> range) {
      DIALS_ASSERT(expr[0] >= 0);
      DIALS_ASSERT(expr[1] > expr[0]);
      DIALS_ASSERT(range[1] > range[0]);
      DIALS_ASSERT(index[1] > index[0]);
      if (groups_.size() > 0) {
        DIALS_ASSERT(expr[0] == groups_[groups_.size() - 1].expr()[1]);
        DIALS_ASSERT(index[0] == groups_[groups_.size() - 1].index()[1]);
      } else {
        DIALS_ASSERT(expr[0] == 0);
      }
      groups_.push_back(Group(index, expr, range));
    }

    /**
     * @returns The number of groups
     */
    std::size_t size() const {
      return groups_.size();
    }

    /**
     * @returns The group
     */
    const Group &operator[](std::size_t index) const {
      DIALS_ASSERT(index < size());
      return groups_[index];
    }

  private:
    std::vector<Group> groups_;
  };

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
      Job(std::size_t index, tiny<int, 2> expr, tiny<int, 2> frames)
          : index_(index), expr_(expr), frames_(frames) {
        DIALS_ASSERT(expr[0] >= 0 && expr[1] > expr[0]);
        DIALS_ASSERT(frames[1] > frames[0]);
      }

      /** @returns The group index */
      std::size_t index() const {
        return index_;
      }

      /** @returns The experiments which this job covers */
      tiny<int, 2> expr() const {
        return expr_;
      }

      /** @returns The frames which this job covers */
      tiny<int, 2> frames() const {
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
      std::size_t index_;
      tiny<int, 2> expr_;
      tiny<int, 2> frames_;
    };

    JobList() {}

    JobList(tiny<int, 2> expr, const af::const_ref<tiny<int, 2> > &jobs) {
      DIALS_ASSERT(expr[1] > expr[0]);
      DIALS_ASSERT(jobs.size() > 0);
      DIALS_ASSERT(jobs[0][1] > jobs[0][0]);
      jobs_.push_back(Job(0, expr, jobs[0]));
      for (std::size_t i = 1; i < jobs.size(); ++i) {
        DIALS_ASSERT(jobs[i][1] > jobs[i][0]);
        DIALS_ASSERT(jobs[i][0] > jobs[i - 1][0]);
        DIALS_ASSERT(jobs[i][1] > jobs[i - 1][1]);
        DIALS_ASSERT(jobs[i][0] <= jobs[i - 1][1]);
        jobs_.push_back(Job(0, expr, jobs[i]));
      }
      groups_.add(int2(0, size()), expr, int2(jobs.front()[0], jobs.back()[1]));
    }

    /**
     * Add a new group of jobs covering a range of experiments
     * @param expr The range of experiments
     * @param range The range of frames
     * @param block_size The job block size
     */
    void add(tiny<int, 2> expr, tiny<int, 2> range, int block_size) {
      std::size_t j0 = size();
      add_jobs(groups_.size(), expr, range, block_size);
      std::size_t j1 = size();
      groups_.add(int2(j0, j1), expr, range);
    }

    /**
     * @returns The requested job
     */
    const Job &operator[](std::size_t index) const {
      DIALS_ASSERT(index < jobs_.size());
      return jobs_[index];
    }

    /**
     * @returns The number of jobs
     */
    std::size_t size() const {
      return jobs_.size();
    }

    /**
     * @returns The group list
     */
    const GroupList groups() const {
      return groups_;
    }

  private:
    void add_jobs(std::size_t index,
                  tiny<int, 2> expr,
                  tiny<int, 2> range,
                  int block_size) {
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
          jobs_.push_back(Job(index, expr, tiny<int, 2>(f, f + 1)));
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
          int i2 = indices[i + 2];
          DIALS_ASSERT(i2 > i1);
          jobs_.push_back(Job(index, expr, tiny<int, 2>(i1, i2)));
        }
        DIALS_ASSERT(jobs_.size() > 0);
      }
    }

    std::vector<Job> jobs_;
    GroupList groups_;
  };

  /**
   * Helper class for checking the range of possible jobs for a particular
   * frame.
   */
  class JobRangeLookup {
  public:
    /**
     * Construct the lookup
     */
    JobRangeLookup(const JobList &jobs) {
      const GroupList &groups = jobs.groups();
      DIALS_ASSERT(0 == groups[0].expr()[0]);
      for (std::size_t i = 0; i < groups.size(); ++i) {
        for (std::size_t j = groups[i].expr()[0]; j < groups[i].expr()[1]; ++j) {
          group_.push_back(i);
        }
      }
      DIALS_ASSERT(group_.size() == groups[groups.size() - 1].expr()[1]);
      offset_.push_back(0);
      for (std::size_t i = 0; i < groups.size(); ++i) {
        tiny<int, 2> f = groups[i].frames();
        DIALS_ASSERT(f[1] > f[0]);
        frame0_.push_back(f[0]);
        offset_.push_back(offset_.back() + f[1] - f[0]);
      }
      DIALS_ASSERT(offset_.back() > 0);
      lookup0_.resize(offset_.back());
      lookup1_.resize(offset_.back());
      for (std::size_t i = 0; i < groups.size(); ++i) {
        std::size_t job0 = groups[i].index()[0];
        std::size_t job1 = groups[i].index()[1];
        DIALS_ASSERT(job1 > job0 && job1 <= jobs.size());
        std::size_t off0 = offset_[i];
        std::size_t off1 = offset_[i + 1];
        DIALS_ASSERT(off1 > off0 && off1 <= lookup0_.size());
        int frame0 = groups[i].frames()[0];
        int frame1 = groups[i].frames()[1];
        DIALS_ASSERT(frame1 > frame0);
        DIALS_ASSERT(frame1 - frame0 == off1 - off0);
        int frame = frame0;
        for (std::size_t i = job0; i < job1; ++i) {
          tiny<int, 2> b = jobs[i].frames();
          DIALS_ASSERT(frame >= b[0]);
          for (; frame < b[1]; ++frame) {
            lookup0_[off0 + frame - frame0] = i;
          }
        }
        DIALS_ASSERT(frame == frame1);
        for (std::size_t i = job1; i > job0; --i) {
          tiny<int, 2> b = jobs[i - 1].frames();
          DIALS_ASSERT(frame <= b[1]);
          for (; frame > b[0]; --frame) {
            lookup1_[off0 + frame - frame0 - 1] = i - 1;
          }
        }
        DIALS_ASSERT(frame == frame0);
        for (std::size_t i = off0 + 1; i < off1; ++i) {
          DIALS_ASSERT(lookup0_[i] >= lookup0_[i - 1]);
          DIALS_ASSERT(lookup1_[i] >= lookup1_[i - 1]);
        }
      }
    }

    /**
     * Get the first job index
     */
    std::size_t first(std::size_t id, int frame) const {
      DIALS_ASSERT(id < group_.size());
      std::size_t group = group_[id];
      DIALS_ASSERT(group < offset_.size() - 1);
      std::size_t offset = offset_[group];
      int frame0 = frame0_[group];
      DIALS_ASSERT(frame >= frame0);
      DIALS_ASSERT(frame < frame0 + (int)offset_[group + 1]);
      std::size_t index = offset + (std::size_t)(frame - frame0);
      DIALS_ASSERT(index < lookup0_.size());
      return lookup0_[index];
    }

    /**
     * Get the second job index
     */
    std::size_t last(std::size_t id, int frame) const {
      DIALS_ASSERT(id < group_.size());
      std::size_t group = group_[id];
      DIALS_ASSERT(group < offset_.size() - 1);
      std::size_t offset = offset_[group];
      int frame0 = frame0_[group];
      DIALS_ASSERT(frame >= frame0);
      DIALS_ASSERT(frame < frame0 + (int)offset_[group + 1]);
      std::size_t index = offset + (std::size_t)(frame - frame0);
      DIALS_ASSERT(index < lookup1_.size());
      return lookup1_[index];
    }

  private:
    std::vector<std::size_t> lookup0_;
    std::vector<std::size_t> lookup1_;
    std::vector<std::size_t> offset_;
    std::vector<std::size_t> group_;
    std::vector<int> frame0_;
  };

  /**
   * A class to managing reflection lookup indices
   */
  class ReflectionLookup {
  public:
    ReflectionLookup(const af::const_ref<int> &id,
                     const af::const_ref<std::size_t> &flags,
                     const af::const_ref<int6> &bbox,
                     const JobList &jobs)
        : jobs_(jobs) {
      DIALS_ASSERT(jobs_.size() > 0);

      typedef std::vector<std::size_t> job_list_type;

      // Check all the reflections are in range
      for (std::size_t i = 0; i < bbox.size(); ++i) {
        DIALS_ASSERT(id[i] >= 0);
        DIALS_ASSERT(bbox[i][1] > bbox[i][0]);
        DIALS_ASSERT(bbox[i][3] > bbox[i][2]);
        DIALS_ASSERT(bbox[i][5] > bbox[i][4]);
      }

      // Compute the job range lookup table
      JobRangeLookup lookup(jobs);

      // Get which reflections to process in which job and task
      std::vector<job_list_type> indices(jobs_.size());
      for (std::size_t index = 0; index < bbox.size(); ++index) {
        std::size_t eid = id[index];
        int z0 = bbox[index][4];
        int z1 = bbox[index][5];
        const std::size_t &f = flags[index];
        if (!(f & af::DontIntegrate)) {
          std::size_t j0 = lookup.first(eid, z0);
          std::size_t j1 = lookup.last(eid, z1 - 1);
          DIALS_ASSERT(j0 < jobs_.size());
          DIALS_ASSERT(j1 < jobs_.size());
          DIALS_ASSERT(j1 >= j0);
          DIALS_ASSERT(z0 >= jobs_[j0].frames()[0]);
          DIALS_ASSERT(z1 <= jobs_[j1].frames()[1]);
          std::size_t jmin = 0;
          double dmin = 0;
          bool inside = false;
          for (std::size_t j = j0; j <= j1; ++j) {
            int jz0 = jobs_[j].frames()[0];
            int jz1 = jobs_[j].frames()[1];
            if (z0 >= jz0 && z1 <= jz1) {
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
          int jz0 = jobs_[jmin].frames()[0];
          int jz1 = jobs_[jmin].frames()[1];
          DIALS_ASSERT(inside == true);
          DIALS_ASSERT(z0 >= jz0 && z1 <= jz1);
          indices[jmin].push_back(index);
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
      std::size_t k = 0;
      for (std::size_t i = 0; i < indices.size(); ++i) {
        const job_list_type &ind = indices[i];
        for (std::size_t j = 0; j < ind.size(); ++j, ++k) {
          DIALS_ASSERT(k < indices_.size());
          indices_[k] = ind[j];
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
    const JobList::Job &job(std::size_t index) const {
      DIALS_ASSERT(index < jobs_.size());
      return jobs_[index];
    }

    /**
     * Get the indices for each job
     */
    af::const_ref<std::size_t> indices(std::size_t index) const {
      DIALS_ASSERT(index < offset_.size() - 1);
      std::size_t i0 = offset_[index];
      std::size_t i1 = offset_[index + 1];
      DIALS_ASSERT(i1 >= i0);
      std::size_t off = i0;
      std::size_t num = i1 - i0;
      DIALS_ASSERT(off + num <= indices_.size());
      return af::const_ref<std::size_t>(&indices_[off], num);
    }

    JobList jobs_;
    af::shared<std::size_t> offset_;
    af::shared<std::size_t> indices_;
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
    ReflectionManager(const JobList &jobs, af::reflection_table data)
        : lookup_(init(jobs, data)), data_(data), finished_(lookup_.size(), false) {
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
    const JobList::Job &job(std::size_t index) const {
      return lookup_.job(index);
    }

    /**
     * @returns The number of reflections in a job
     */
    std::size_t num_reflections(std::size_t index) const {
      DIALS_ASSERT(index < finished_.size());
      return lookup_.indices(index).size();
    }

    /**
     * @returns The reflections for a particular block.
     */
    af::reflection_table split(std::size_t index) {
      using namespace af::boost_python::flex_table_suite;

      // Check the input
      DIALS_ASSERT(index < finished_.size());
      af::const_ref<std::size_t> ind = lookup_.indices(index);

      // Extract the reflection table
      af::reflection_table result = select_rows_index(data_, ind);

      // Make sure that the experiment ids start at zero
      /* tiny<int,2> expr = job(index).expr(); */
      /* if (expr[0] != 0) { */
      /*   DIALS_ASSERT(expr[0] > 0); */
      /*   DIALS_ASSERT(expr[1] > expr[0]); */
      /*   af::ref<std::size_t> id = result["id"]; */
      /*   for (std::size_t i = 0; i < id.size(); ++i) { */
      /*     DIALS_ASSERT(id[i] >= expr[0]); */
      /*     DIALS_ASSERT(id[i] < expr[1]); */
      /*     id[i] -= expr[0]; */
      /*   } */
      /* } */

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
      DIALS_ASSERT(ind.size() == result.size());

      // Rejig the experiment ids again
      /* tiny<int,2> expr = job(index).expr(); */
      /* DIALS_ASSERT(expr[0] >= 0); */
      /* DIALS_ASSERT(expr[1] > expr[0]); */
      /* std::size_t num_expr = expr[1] - expr[0]; */
      /* af::ref<std::size_t> id = result["id"]; */
      /* for (std::size_t i = 0; i < id.size(); ++i) { */
      /*   DIALS_ASSERT(id[i] < num_expr); */
      /*   id[i] += expr[0]; */
      /* } */

      // Set the result
      set_selected_rows_index(data_, ind, result);

      // Set finished flag
      finished_[index] = true;
    }

  private:
    /**
     * Initialise the indexer
     */
    ReflectionLookup init(const JobList &jobs, af::reflection_table data) {
      DIALS_ASSERT(data.is_consistent());
      DIALS_ASSERT(data.size() > 0);
      DIALS_ASSERT(data.contains("id"));
      DIALS_ASSERT(data.contains("flags"));
      DIALS_ASSERT(data.contains("bbox"));
      DIALS_ASSERT(jobs.size() > 0);
      return ReflectionLookup(data["id"], data["flags"], data["bbox"], jobs);
    }

    ReflectionLookup lookup_;
    af::reflection_table data_;
    af::shared<bool> finished_;
  };

}}  // namespace dials::algorithms

#endif  // DIALS_ALGORITHMS_INTEGRATION_INTEGRATOR_H
