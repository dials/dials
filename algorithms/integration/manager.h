/*
 * manager.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_INTEGRATION_MANAGER_H
#define DIALS_ALGORITHMS_INTEGRATION_MANAGER_H

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

  /**
   * A class to managing reflection lookup indices
   */
  class ReflectionLookup2 {
  public:
    ReflectionLookup2(const af::const_ref<int> &id,
                      const af::const_ref<int6> &bbox,
                      const int2 &frames) {
      DIALS_ASSERT(frames[1] > frames[0]);

      // Make a list of the frames
      for (int i = frames[0]; i < frames[1]; ++i) {
        frames_.push_back(int2(i, i + 1));
      }

      // Check all the reflection bboxes are valid
      for (std::size_t i = 0; i < bbox.size(); ++i) {
        DIALS_ASSERT(id[i] >= 0);
        DIALS_ASSERT(bbox[i][1] > bbox[i][0]);
        DIALS_ASSERT(bbox[i][3] > bbox[i][2]);
        DIALS_ASSERT(bbox[i][5] > bbox[i][4]);
        DIALS_ASSERT(bbox[i][4] >= frames[0]);
        DIALS_ASSERT(bbox[i][5] <= frames[1]);
      }

      // Count the number of reflections on each image
      std::vector<int> count(frames_.size(), 0);
      for (std::size_t i = 0; i < bbox.size(); ++i) {
        for (int j = bbox[i][4]; j < bbox[i][5]; ++j) {
          int k = j - frames[0];
          DIALS_ASSERT(k >= 0);
          DIALS_ASSERT(k < count.size());
          count[k]++;
        }
      }

      // Compute offsests
      offset_.push_back(0);
      std::partial_sum(count.begin(), count.end(), std::back_inserter(offset_));
      DIALS_ASSERT(offset_.size() == count.size() + 1);

      // Compute indices
      indices_.resize(offset_.back());
      count.assign(count.size(), 0);
      for (std::size_t i = 0; i < bbox.size(); ++i) {
        for (int j = bbox[i][4]; j < bbox[i][5]; ++j) {
          int k = j - frames[0];
          DIALS_ASSERT(k >= 0);
          DIALS_ASSERT(k < count.size());
          std::size_t l = offset_[k] + count[k]++;
          DIALS_ASSERT(l < offset_[k + 1]);
          DIALS_ASSERT(l < indices_.size());
          indices_[l] = i;
        }
      }
      DIALS_ASSERT(indices_.size() == bbox.size());

      // Check we didn't make any mistakes
      for (std::size_t i = 0; i < count.size(); ++i) {
        DIALS_ASSERT(offset_[i] + count[i] == offset_[i + 1]);
      }
    }

    /**
     * @returns The number of tasks
     */
    std::size_t size() const {
      return frames_.size();
    }

    /**
     * @returns The block indices
     */
    int2 frames(std::size_t index) const {
      DIALS_ASSERT(index < frames_.size());
      return frames_[index];
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

    af::shared<int2> frames_;
    af::shared<std::size_t> offset_;
    af::shared<std::size_t> indices_;
  };

  class ReflectionManagerPerImage {
  public:
    ReflectionManagerPerImage(int2 frames, af::reflection_table data)
        : lookup_(init(frames, data)), data_(data), finished_(lookup_.size(), false) {
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
    int2 frames(std::size_t index) const {
      return lookup_.frames(index);
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
      DIALS_ASSERT(index < finished_.size());
      return select_rows_index(data_, lookup_.indices(index));
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

      // Set the result
      set_selected_rows_index(data_, ind, result);

      // Set finished flag
      finished_[index] = true;
    }

  private:
    /**
     * Initialise the indexer
     */
    ReflectionLookup2 init(int2 frames, af::reflection_table data) {
      DIALS_ASSERT(data.is_consistent());
      DIALS_ASSERT(data.size() > 0);
      DIALS_ASSERT(data.contains("id"));
      DIALS_ASSERT(data.contains("bbox"));
      DIALS_ASSERT(frames[1] > frames[0]);
      return ReflectionLookup2(data["id"], data["bbox"], frames);
    }

    ReflectionLookup2 lookup_;
    af::reflection_table data_;
    af::shared<bool> finished_;
  };

}}  // namespace dials::algorithms

#endif  // DIALS_ALGORITHMS_INTEGRATION_MANAGER_H
