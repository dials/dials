#ifndef DIALS_ARRAY_FAMILY_BOOST_PYTHON_REFLECTION_TABLE_SUITE_H
#define DIALS_ARRAY_FAMILY_BOOST_PYTHON_REFLECTION_TABLE_SUITE_H

#include <dxtbx/array_family/flex_table_suite.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/error.h>
#include <dxtbx/model/experiment.h>
#include <dxtbx/model/experiment_list.h>

namespace dials { namespace af { namespace boost_python {
  namespace reflection_table_suite {

    namespace flex_table_suite = dxtbx::af::flex_table_suite;

    /**
     * Select a number of rows from the table via an index array
     * and copy across experiment identifiers
     * @param self The current table
     * @param index The index array
     * @returns The new table with the requested rows
     */
    template <typename T>
    T select_rows_index_base(const T &self,
                        const scitbx::af::const_ref<std::size_t> &index,
                        const bool preserve_ids) {
      T new_table = flex_table_suite::select_rows_index(self, index);

      // Get the id column (if it exists) and make a set of unique values
      if (preserve_ids && self.contains("id")) {
        af::shared<int> col = new_table["id"];
        std::set<int> new_ids(col.begin(), col.end());

        // Copy across identifiers for ids in new table
        typedef typename T::experiment_map_type::const_iterator const_iterator;
        for (std::set<int>::iterator i = new_ids.begin(); i != new_ids.end(); ++i) {
          const_iterator found = self.experiment_identifiers()->find(*i);
          if (found != self.experiment_identifiers()->end()) {
            (*new_table.experiment_identifiers())[found->first] = found->second;
          }
        }
      }
      return new_table;
    }

    template <typename T>
    T select_rows_index(const T &self,
                        const scitbx::af::const_ref<std::size_t> &index) {
      return select_rows_index_base(self, index, true);
    }
    template <typename T>
    T select_rows_index_fast(const T &self,
                        const scitbx::af::const_ref<std::size_t> &index) {
      return select_rows_index_base(self, index, false);
    }
                             

    /**
     * Select a number of rows from the table via an index array
     * @param self The current table
     * @param flags The flag array
     * @returns The new table with the requested rows
     */
    template <typename T>
    T select_rows_flags(const T &self, const af::const_ref<bool> &flags) {
      DIALS_ASSERT(self.nrows() == flags.size());
      af::shared<std::size_t> index;
      index.reserve(flags.size());
      for (std::size_t i = 0; i < flags.size(); ++i) {
        if (flags[i]) index.push_back(i);
      }
      return select_rows_index(self, index.const_ref());
    }

    template <typename T>
    T select_rows_flags_fast(const T &self, const af::const_ref<bool> &flags) {
      DIALS_ASSERT(self.nrows() == flags.size());
      af::shared<std::size_t> index;
      index.reserve(flags.size());
      for (std::size_t i = 0; i < flags.size(); ++i) {
        if (flags[i]) index.push_back(i);
      }
      return select_rows_index_fast(self, index.const_ref());
    }


    /**
     * Adds missing experiment identifiers to self that are in other
     * @param self The current table
     * @param other The table to add experiment identifiers from
     */
    template <typename T>
    void extend_identifiers(T &self, const T &other) {
      typedef typename T::experiment_map_type::const_iterator const_iterator;
      typedef typename T::experiment_map_type::iterator iterator;
      for (const_iterator it = other.experiment_identifiers()->begin();
           it != other.experiment_identifiers()->end();
           ++it) {
        iterator found = self.experiment_identifiers()->find(it->first);
        if (found == self.experiment_identifiers()->end()) {
          (*self.experiment_identifiers())[it->first] = it->second;
        } else if (it->second != found->second) {
          throw DIALS_ERROR("Experiment identifiers do not match");
        }
      }
    }

    /**
     * Extend the table with column data from another table. This will add
     * all the columns from the other table onto the end of the current table,
     * and add any missing experiment identifiers.
     * @param self The current table
     * @param other The other table
     */
    template <typename T>
    void extend(T &self, const T &other) {
      flex_table_suite::extend(self, other);
      extend_identifiers(self, other);
    }

    /**
     * Select a number of rows from the table using an experiment identifier
     * @param self The current table
     * @param expt The experiment containing the identifier for selection
     * @returns A new table containing rows that match the experiment identifier
     */
    template <typename T>
    T select_using_experiment(T &self, dxtbx::model::Experiment expt) {
      typedef typename T::experiment_map_type::const_iterator const_iterator;

      std::string identifier = expt.get_identifier();
      int id_value = -1;
      for (const_iterator it = self.experiment_identifiers()->begin();
           it != self.experiment_identifiers()->end();
           ++it) {
        if (identifier == it->second) {
          id_value = it->first;
          break;
        }
      }

      T result;
      if (self.contains("id") && id_value != -1) {
        af::shared<int> col1 = self["id"];
        af::shared<std::size_t> sel;
        for (int i = 0; i < col1.size(); ++i) {
          if (col1[i] == id_value) {
            sel.push_back(i);
          }
        }
        af::const_ref<std::size_t> idx = sel.const_ref();

        result = select_rows_index(self, idx);
      }

      return result;
    }

    /**
     * Select a number of rows from the table using experiment identifiers
     * @param self The current table
     * @param expts The experiments containing the identifiers used for selection
     * @returns A new table containing rows that match the experiment identifiers
     */
    template <typename T>
    T select_using_experiments(T &self, dxtbx::model::ExperimentList expts) {
      typedef typename T::experiment_map_type::const_iterator const_iterator;
      typedef dxtbx::model::ExperimentList::shared_type::const_iterator
        expt_const_iterator;
      T result;
      for (expt_const_iterator expt = expts.begin(); expt != expts.end(); ++expt) {
        std::string identifier = expt->get_identifier();
        int id_value = -1;
        for (const_iterator it = self.experiment_identifiers()->begin();
             it != self.experiment_identifiers()->end();
             ++it) {
          if (identifier == it->second) {
            id_value = it->first;
            break;
          }
        }

        if (self.contains("id") && id_value != -1) {
          af::shared<int> col1 = self["id"];
          af::shared<std::size_t> sel;
          for (int i = 0; i < col1.size(); ++i) {
            if (col1[i] == id_value) {
              sel.push_back(i);
            }
          }
          af::const_ref<std::size_t> idx = sel.const_ref();

          T sel_refl = select_rows_index(self, idx);
          extend(result, sel_refl);
        }
      }

      return result;
    }

    /**
     * Update the table with column data from another table. New columns are added
     * to the table and existing columns are over-written by columns from the
     * other table. Before updating, missing experiment identifiers from other are
     * added to self.
     * @param self The current table
     * @param other The other table
     */
    template <typename T>
    void update(T &self, const T &other) {
      extend_identifiers(self, other);
      flex_table_suite::update(self, other);
    }

    template <typename T>
    T deepcopy(const T &self, boost::python::dict obj) {
      T new_table = flex_table_suite::deepcopy(self, obj);

      typedef typename T::experiment_map_type::const_iterator const_iterator;
      for (const_iterator it = self.experiment_identifiers()->begin();
           it != self.experiment_identifiers()->end();
           ++it) {
        (*new_table.experiment_identifiers())[it->first] = it->second;
      }
      return new_table;
    }

    /**
     * Select a number of columns from the table via a key array
     * @param self The current table
     * @param keys The key array
     * @returns The new table with the requested columns
     */
    template <typename T>
    T select_cols_keys(const T &self, const af::const_ref<std::string> &keys) {
      T result = flex_table_suite::select_cols_keys<T>(self, keys);
      extend_identifiers(result, self);
      return result;
    }

    /**
     * Select a number of columns from the table via an key array
     * @param self The current table
     * @param keys The key array
     * @returns The new table with the requested columns
     */
    template <typename T>
    T select_cols_tuple(const T &self, boost::python::tuple keys) {
      T result = flex_table_suite::select_cols_tuple<T>(self, keys);
      extend_identifiers(result, self);
      return result;
    }

    /**
     * Get a slice of the table and return a new table
     * @param self The current table
     * @param slice The slice
     * @returns A new table with the chosen elements
     */
    template <typename T>
    T getitem_slice(const T &self, boost::python::slice s) {
      typedef typename T::const_iterator iterator;
      scitbx::boost_python::adapted_slice as(s, self.nrows());
      T result(as.size);
      for (iterator it = self.begin(); it != self.end(); ++it) {
        dxtbx::af::flex_table_suite::copy_to_slice_visitor<T> visitor(
          result, it->first, as);
        it->second.apply_visitor(visitor);
      }
      if (self.contains("id")) {
        /* note some tables contain id values of -1 for unindexed reflections
        but the identifiers map only allows keys of type size_t
        */
        af::shared<int> col = result["id"];
        std::set<int> new_ids(col.begin(), col.end());
        typedef typename T::experiment_map_type::iterator iterator;
        for (std::set<int>::iterator i = new_ids.begin(); i != new_ids.end(); ++i) {
          iterator found = self.experiment_identifiers()->find(*i);
          if (found != self.experiment_identifiers()->end()) {
            (*result.experiment_identifiers())[found->first] = found->second;
          }
        }
      }
      return result;
    }

}}}}  // namespace dials::af::boost_python::reflection_table_suite

#endif  // DIALS_ARRAY_FAMILY_BOOST_PYTHON_REFLECTION_TABLE_SUITE_H
