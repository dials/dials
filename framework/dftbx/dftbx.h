/* dftbx
 *
 * DIALS Framework Toolbox:
 *
 * Framework code (i.e. not model, not algorithm, not interface) code for DIALS
 * which currently contains:
 *
 * ObservationList template class
 * ReflectionList template class
 *
 */

#ifndef DFTBX_DFTBX_H
#define DFTBX_DFTBX_H

#include <boost/python.hpp>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/flex_types.h>
#include <miller.h>
#include <uctbx.h>
#include <cctype>
#include <map>
#include <algorithm>

namespace dftbx {

  void init_module() { };

  /** miller_index_t
   *
   * Typedef'd to be integer Miller indices: DIALS won't be using other types.
   */

  typedef cctbx::miller::index<int> miller_index_t;

  /** ObservationList<T>
   *
   * If T is big the caller should be responsible for making T a shared pointer
   * and handling issues with duplicate pointers possibly modifying T's
   * elsewhere. If T is small no worries, but the merge() method below will

   *
   */

  template <class T>
    class ObservationList {
  public:

    typedef typename scitbx::af::shared<T> vector_type;

    size_t size() const {
      return observations_.size();
    }

    void push_back(const T & t) {
      observations_.push_back(t);
    }

    void clear() {
      observations_.clear();
    }

    T & operator[](size_t j) {
      return observations_[j];
    }

    const T & operator[](size_t j) const {
      return observations_[j];
    }

    void duplicate(const ObservationList<T> & other) {
      for (size_t j = 0; j < other.size(); j ++) {
        observations_.push_back(other[j]);
      }
    }

    void merge(ObservationList<T> & other) {
      for (size_t j = 0; j < other.size(); j ++) {
        observations_.push_back(other[j]);
      }
      other.clear();
    }

  private:

    vector_type observations_;

  };

  template <class T>
    class ReflectionList {
  public:

    typedef ObservationList<T> observation_list_t;

    typedef std::map<miller_index_t, observation_list_t> map_type;
    typedef typename map_type::iterator map_iterator_t;

    observation_list_t & operator[](miller_index_t hkl) {
      return reflections_[hkl];
    }

    const observation_list_t & operator[](miller_index_t hkl) const {
      return reflections_[hkl];
    }

    void merge(ReflectionList<T> & other) {
      for (map_iterator_t mi = reflections_.begin();
           mi != reflections_.end(); ++mi) {
        mi->second.merge(other[mi->first]);
      }
    }

  private:
    map_type reflections_;

  };

  typedef ObservationList<int> int_ol;
  typedef ReflectionList<int> int_rl;


}

#endif
