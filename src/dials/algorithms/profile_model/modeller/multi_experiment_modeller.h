/*
 * multi_experiment_modeller.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */

#ifndef DIALS_ALGORITHMS_PROFILE_MODEL_MODELLER_MULTI_EXPERIMENT_MODELLER_H
#define DIALS_ALGORITHMS_PROFILE_MODEL_MODELLER_MULTI_EXPERIMENT_MODELLER_H

#include <vector>
#include <numeric>
#include <dials/array_family/reflection_table.h>
#include <dials/array_family/boost_python/flex_table_suite.h>
#include <dials/algorithms/profile_model/modeller/modeller_interface.h>

namespace dials { namespace algorithms {

  /**
   * A class to compute reference profiles for each experiment.
   */
  class MultiExpProfileModeller {
  public:
    typedef ProfileModellerIface::pointer modeller_pointer;

    /**
     * Initialize the multi experiment profile modeller
     */
    MultiExpProfileModeller() {}

    /**
     * Add a profile modeller
     * @param modeller The profile modeller
     */
    void add(modeller_pointer modeller) {
      DIALS_ASSERT(modeller != NULL);
      modellers_.push_back(modeller);
    }

    /**
     * Get the modeller for a particular index
     * @param index The index
     * @return The modeller
     */
    modeller_pointer operator[](std::size_t index) const {
      DIALS_ASSERT(index < modellers_.size());
      return modellers_[index];
    }

    /**
     * Model the reflections
     * @param reflections The reflection table
     */
    void model(af::reflection_table reflections) {
      using af::boost_python::flex_table_suite::select_rows_index;
      using af::boost_python::flex_table_suite::set_selected_rows_index;

      // Check some stuff
      DIALS_ASSERT(size() > 0);
      DIALS_ASSERT(reflections.size() > 0);
      DIALS_ASSERT(reflections.contains("id"));

      // Get the experiment id
      af::const_ref<int> id = reflections["id"];

      // Compute the number of reflections for each experiment
      std::vector<std::size_t> num1(size(), 0);
      for (std::size_t i = 0; i < id.size(); ++i) {
        DIALS_ASSERT(id[i] >= 0);
        DIALS_ASSERT(id[i] < num1.size());
        num1[id[i]]++;
      }

      // Compute the offset array
      std::vector<std::size_t> offset(1, 0);
      std::partial_sum(num1.begin(), num1.end(), std::back_inserter(offset));
      DIALS_ASSERT(offset.size() == num1.size() + 1);
      DIALS_ASSERT(offset.back() == id.size());

      // Compute the indices
      std::vector<std::size_t> indices(id.size());
      std::vector<std::size_t> num2(num1.size(), 0);
      for (std::size_t i = 0; i < id.size(); ++i) {
        DIALS_ASSERT(id[i] < num1.size());
        std::size_t o1 = offset[id[i]];
        std::size_t o2 = offset[id[i] + 1];
        std::size_t n1 = num1[id[i]];
        std::size_t n2 = num2[id[i]];
        std::size_t j = o1 + n2;
        DIALS_ASSERT(j < o2);
        DIALS_ASSERT(n2 < n1);
        DIALS_ASSERT(j < indices.size());
        indices[j] = i;
        num2[id[i]]++;
      }

      // Check we've assigned everything
      for (std::size_t i = 0; i < num1.size(); ++i) {
        DIALS_ASSERT(num1[i] == num2[i]);
      }

      // Process all the reflections
      for (std::size_t i = 0; i < modellers_.size(); ++i) {
        // Get the indices
        std::size_t o1 = offset[i];
        std::size_t o2 = offset[i + 1];
        DIALS_ASSERT(o2 <= indices.size());
        DIALS_ASSERT(o1 <= o2);
        std::size_t n = o2 - o1;

        // If there are reflections for this experiment do some modelling
        if (n > 0) {
          af::const_ref<std::size_t> ind(&indices[o1], n);
          af::reflection_table subset = select_rows_index(reflections, ind);
          DIALS_ASSERT(modellers_[i] != NULL);
          modellers_[i]->model(subset);
          set_selected_rows_index(reflections, ind, subset);
        }
      }
    }

    /**
     * Do the profile fitting
     * @param reflections
     */
    af::shared<bool> fit(af::reflection_table reflections) const {
      using af::boost_python::flex_table_suite::select_rows_index;
      using af::boost_python::flex_table_suite::set_selected_rows_index;

      // Check some stuff
      DIALS_ASSERT(size() > 0);
      DIALS_ASSERT(reflections.size() > 0);
      DIALS_ASSERT(reflections.contains("id"));

      // Get the experiment id
      af::const_ref<int> id = reflections["id"];

      // Compute the number of reflections for each experiment
      std::vector<std::size_t> num1(size(), 0);
      for (std::size_t i = 0; i < id.size(); ++i) {
        DIALS_ASSERT(id[i] >= 0);
        DIALS_ASSERT(id[i] < num1.size());
        num1[id[i]]++;
      }

      // Compute the offset array
      std::vector<std::size_t> offset(1, 0);
      std::partial_sum(num1.begin(), num1.end(), std::back_inserter(offset));
      DIALS_ASSERT(offset.size() == num1.size() + 1);
      DIALS_ASSERT(offset.back() == id.size());

      // Compute the indices
      std::vector<std::size_t> indices(id.size());
      std::vector<std::size_t> num2(num1.size(), 0);
      for (std::size_t i = 0; i < id.size(); ++i) {
        DIALS_ASSERT(id[i] < num1.size());
        std::size_t o1 = offset[id[i]];
        std::size_t o2 = offset[id[i] + 1];
        std::size_t n1 = num1[id[i]];
        std::size_t n2 = num2[id[i]];
        std::size_t j = o1 + n2;
        DIALS_ASSERT(j < o2);
        DIALS_ASSERT(n2 < n1);
        indices[j] = i;
        num2[id[i]]++;
      }

      // Check we've assigned everything
      for (std::size_t i = 0; i < num1.size(); ++i) {
        DIALS_ASSERT(num1[i] == num2[i]);
      }

      // Process all the reflections
      af::shared<bool> success(reflections.size(), false);
      for (std::size_t i = 0; i < modellers_.size(); ++i) {
        // Get the indices
        std::size_t o1 = offset[i];
        std::size_t o2 = offset[i + 1];
        DIALS_ASSERT(o2 <= indices.size());
        DIALS_ASSERT(o1 <= o2);
        std::size_t n = o2 - o1;

        // If there are any reflections, do the fitting
        if (n > 0) {
          af::const_ref<std::size_t> ind(&indices[o1], n);
          af::reflection_table subset = select_rows_index(reflections, ind);
          af::shared<bool> subset_success = modellers_[i]->fit(subset);
          set_selected_rows_index(reflections, ind, subset);
          for (std::size_t j = 0; j < ind.size(); ++j) {
            success[ind[j]] = subset_success[j];
          }
        }
      }

      // Return success
      return success;
    }

    /**
     * Do the profile validation
     * @param reflections
     */
    void validate(af::reflection_table reflections) const {
      using af::boost_python::flex_table_suite::select_rows_index;
      using af::boost_python::flex_table_suite::set_selected_rows_index;

      // Check some stuff
      DIALS_ASSERT(size() > 0);
      DIALS_ASSERT(reflections.size() > 0);
      DIALS_ASSERT(reflections.contains("id"));

      // Get the experiment id
      af::const_ref<int> id = reflections["id"];

      // Compute the number of reflections for each experiment
      std::vector<std::size_t> num1(size(), 0);
      for (std::size_t i = 0; i < id.size(); ++i) {
        DIALS_ASSERT(id[i] >= 0);
        DIALS_ASSERT(id[i] < num1.size());
        num1[id[i]]++;
      }

      // Compute the offset array
      std::vector<std::size_t> offset(1, 0);
      std::partial_sum(num1.begin(), num1.end(), std::back_inserter(offset));
      DIALS_ASSERT(offset.size() == num1.size() + 1);
      DIALS_ASSERT(offset.back() == id.size());

      // Compute the indices
      std::vector<std::size_t> indices(id.size());
      std::vector<std::size_t> num2(num1.size(), 0);
      for (std::size_t i = 0; i < id.size(); ++i) {
        DIALS_ASSERT(id[i] < num1.size());
        std::size_t o1 = offset[id[i]];
        std::size_t o2 = offset[id[i] + 1];
        std::size_t n1 = num1[id[i]];
        std::size_t n2 = num2[id[i]];
        std::size_t j = o1 + n2;
        DIALS_ASSERT(j < o2);
        DIALS_ASSERT(n2 < n1);
        indices[j] = i;
        num2[id[i]]++;
      }

      // Check we've assigned everything
      for (std::size_t i = 0; i < num1.size(); ++i) {
        DIALS_ASSERT(num1[i] == num2[i]);
      }

      // Process all the reflections
      for (std::size_t i = 0; i < modellers_.size(); ++i) {
        // Get the indices
        std::size_t o1 = offset[i];
        std::size_t o2 = offset[i + 1];
        DIALS_ASSERT(o2 <= indices.size());
        DIALS_ASSERT(o1 < o2);
        std::size_t n = o2 - o1;

        // The indices
        af::const_ref<std::size_t> ind(&indices[o1], n);

        // Get the reflections
        af::reflection_table subset = select_rows_index(reflections, ind);

        // Do the fitting
        modellers_[i]->validate(subset);

        // Set any results
        set_selected_rows_index(reflections, ind, subset);
      }
    }

    /**
     * Add the profiles from another modeller
     * @param other The other modeller
     */
    void accumulate(const MultiExpProfileModeller &other) {
      DIALS_ASSERT(size() == other.size());
      for (std::size_t i = 0; i < size(); ++i) {
        modellers_[i]->accumulate(other.modellers_[i]);
      }
    }

    /**
     * Finalize the profiles
     */
    void finalize() {
      for (std::size_t i = 0; i < modellers_.size(); ++i) {
        modellers_[i]->finalize();
      }
    }

    /**
     * @return Is the model finalized
     */
    bool finalized() const {
      bool result = false;
      for (std::size_t i = 0; i < modellers_.size(); ++i) {
        if (modellers_[i]->finalized()) {
          result = true;
          break;
        }
      }
      return result;
    }

    /**
     * @return the number of modllers.
     */
    std::size_t size() const {
      return modellers_.size();
    }

    /**
     * Do a deep copy
     */
    MultiExpProfileModeller copy() const {
      MultiExpProfileModeller result;
      for (std::size_t i = 0; i < modellers_.size(); ++i) {
        result.modellers_.push_back(modellers_[i]->copy());
      }
      return result;
    }

  private:
    std::vector<modeller_pointer> modellers_;
  };

}}  // namespace dials::algorithms

#endif  // DIALS_ALGORITHMS_PROFILE_MODEL_MODELLER_MULTI_EXPERIMENT_MODELLER_H
