/*
 * reflection_table.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ARRAY_FAMILY_REFLECTION_TABLE_H
#define DIALS_ARRAY_FAMILY_REFLECTION_TABLE_H

#include <dials/array_family/flex_table.h>
#include <dials/model/data/shoebox.h>
#include <scitbx/array_family/tiny_types.h>
#include <scitbx/mat3.h>
#include <scitbx/vec3.h>
#include <scitbx/vec2.h>
#include <cctbx/miller.h>

namespace dials { namespace af {

  using scitbx::vec2;
  using scitbx::vec3;
  using scitbx::mat3;
  using scitbx::af::int6;
  using model::Shoebox;


  template <> inline
  vec2<double> init_zero< vec2<double> >() {
    return vec2<double>(0.0,0.0);
  }

  template <> inline
  vec3<double> init_zero< vec3<double> >() {
    return vec3<double>(0.0,0.0,0.0);
  }

  template <> inline
  mat3<double> init_zero< mat3<double> >() {
    return mat3<double>(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
  }

  template <> inline
  int6 init_zero<int6>() {
    return int6(0,0,0,0,0,0);
  }

  template <> inline
  cctbx::miller::index<> init_zero< cctbx::miller::index<> >() {
    return cctbx::miller::index<>(0,0,0);
  }


  typedef flex_type_generator<
    bool,
    int,
    std::size_t,
    double,
    std::string,
    vec2<double>,
    vec3<double>,
    mat3<double>,
    int6,
    cctbx::miller::index<>,
    Shoebox<>
  > reflection_table_type_generator;

  typedef typename reflection_table_type_generator::type reflection_table_types;
  typedef flex_table<reflection_table_types> reflection_table;

  enum Flags {

    // Predicted/Observed
    Predicted        = (1 << 0),
    Observed         = (1 << 1),

    // Use in indexing/refinement
    Indexed          = (1 << 2),
    UsedInRefinement = (1 << 3),
    Strong           = (1 << 5),

    // Role in integration
    ReferenceSpot    = (1 << 6),
    DontIntegrate    = (1 << 7),

    // Integated
    IntegratedSum    = (1 << 8),
    IntegratedPrf    = (1 << 9),
    Integrated       = IntegratedSum | IntegratedPrf,

    // Bad shoebox
    Overloaded       = (1 << 10),
    OverlappedBg     = (1 << 11),
    OverlappedFg     = (1 << 12),
    InPowderRing     = (1 << 13),
    ForegroundIncludesBadPixels = (1 << 14),
    BackgroundIncludesBadPixels = (1 << 15),
    IncludesBadPixels = ForegroundIncludesBadPixels | BackgroundIncludesBadPixels,
    BadShoebox       = Overloaded | OverlappedBg | OverlappedFg | InPowderRing | IncludesBadPixels,

    // Bad spot
    BadSpot = BadShoebox,

    // Profile Modelling
    UsedInModelling  = (1 << 16),

    // Centroid outlier
    CentroidOutlier = (1 << 17),

    // Some Error Codes
    FailedDuringBackgroundModelling = (1 << 18),
    FailedDuringSummation = (1 << 19),
    FailedDuringProfileFitting = (1 << 20),

    // Bad reference
    BadReference = (1 << 21),
  };

  /**
   * A class to represent a reflection
   */
  class Reflection {
  public:

    typedef typename reflection_table_type_generator::data_type data_type;
    typedef std::map<std::string, data_type> map_type;

    typedef typename map_type::key_type key_type;
    typedef typename map_type::mapped_type mapped_type;
    typedef typename map_type::value_type map_value_type;
    typedef typename map_type::iterator iterator;
    typedef typename map_type::const_iterator const_iterator;
    typedef typename map_type::size_type size_type;

    /**
     * Instantiate
     */
    Reflection()
      : data_(boost::make_shared<map_type>()){}

    /**
     * Access a value by key
     * @param key The column name
     * @returns The proxy object to access the value
     */
    const mapped_type& operator[](const key_type &key) const {
      return data_->operator[](key);
    }

    /**
     * Access a value by key
     * @param key The column name
     * @returns The proxy object to access the value
     */
    mapped_type& operator[](const key_type &key) {
      return data_->operator[](key);
    }

    /**
     * Access a value by key
     * @param key The column name
     * @returns The value.
     */
    template <typename T>
    T& get(const key_type &key) {
      iterator it = find(key);
      DIALS_ASSERT(it != end());
      return boost::get<T>(it->second);
    }

    /**
     * Access a value by key
     * @param key The column name
     * @returns The value.
     */
    template <typename T>
    const T& get(const key_type &key) const {
      const_iterator it = find(key);
      DIALS_ASSERT(it != end());
      return boost::get<T>(it->second);
    }

    /** @returns An iterator to the beginning of the column map */
    iterator begin() {
      return data_->begin();
    }

    /** @returns An iterator to the end of the column map */
    iterator end() {
      return data_->end();
    }

    /** @returns A const iterator to the beginning of the column map */
    const_iterator begin() const {
      return data_->begin();
    }

    /** @returns A const iterator to the end of the column map */
    const_iterator end() const {
      return data_->end();
    }

    /** @returns The number of values in the table */
    size_type size() const {
      return data_->size();
    }

    /** @returns Is the table empty */
    bool empty() const {
      return data_->empty();
    }

    /** @returns The number of columns matching the key (0 or 1) */
    size_type count(const key_type &key) const {
      return data_->count(key);
    }

    /**
     * Find a column matching the key
     * @param key The column name
     * @returns An iterator to the column
     */
    iterator find(const key_type &key) {
      return data_->find(key);
    }

    /**
     * Find a column matching the key
     * @param key The column name
     * @returns A const iterator to the column
     */
    const_iterator find(const key_type &key) const {
      return data_->find(key);
    }

    /**
     * Erase a column from the table.
     * @param key The column name
     * @returns The number of columns removed
     */
    size_type erase(const key_type &key) {
      return data_->erase(key);
    }

    /** Clear the table */
    void clear() {
      data_->clear();
    }

    /** @returns Does the table contain the key. */
    bool contains(const key_type &key) const {
      const_iterator it = find(key);
      return it != end();
    }

  protected:

    boost::shared_ptr<map_type> data_;

  };

  namespace detail {

    /**
     * A visitor to convert extract reflection table row to a reflection object
     */
    struct row_to_reflection_visitor :
        public boost::static_visitor<Reflection::data_type> {
      std::size_t n_;
      row_to_reflection_visitor(std::size_t n) : n_(n) {}
      template <typename T>
      Reflection::data_type operator () (T &col) {
        DIALS_ASSERT(n_ < col.size());
        return Reflection::data_type(col[n_]);
      }
    };

    /**
     * A visitor to convert a reflection object to a reflection table row
     */
    struct reflection_to_row_visitor : public boost::static_visitor<void> {
      af::reflection_table table_;
      std::size_t n_;
      Reflection::key_type key_;
      reflection_to_row_visitor(af::reflection_table table, std::size_t n, Reflection::key_type key) :
        table_(table),
        n_(n),
        key_(key) {}
      template <typename T>
      void operator () (const T &item) {
        af::ref<T> col = table_[key_];
        DIALS_ASSERT(n_ < col.size());
        col[n_] = item;
      }
    };

    /**
     * Get a reflection object from the reflection table
     * @param table The reflection table
     * @param index The array index
     * @returns The reflection object
     */
    inline
    Reflection reflection_table_get_reflection(
        af::reflection_table table, std::size_t index) {
      typedef typename af::reflection_table::const_iterator iterator;
      DIALS_ASSERT(index < table.size());
      Reflection result;
      row_to_reflection_visitor visitor(index);
      for (iterator it = table.begin(); it != table.end(); ++it) {
        result[it->first] = it->second.apply_visitor(visitor);
      }
      return result;
    }

    /**
     * Set the reflection object in the reflection table
     * @param table The reflection table
     * @param index The array index
     * @param value The reflection
     */
    inline
    void reflection_table_set_reflection(
        af::reflection_table table, std::size_t index, Reflection value) {
      typedef typename Reflection::const_iterator iterator;
      DIALS_ASSERT(index < table.size());
      for (iterator it = value.begin(); it != value.end(); ++it) {
        reflection_to_row_visitor visitor(table, index, it->first);
        it->second.apply_visitor(visitor);
      }
    }

  }

  /**
   * Convert a reflection table to an array of reflections
   * @param table The reflection table
   * @returns The array of reflections
   */
  inline
  af::shared<Reflection> reflection_table_to_array(af::reflection_table table) {
    af::shared<Reflection> result;
    result.reserve(table.size());
    for (std::size_t i = 0; i < table.size(); ++i) {
      result.push_back(detail::reflection_table_get_reflection(table, i));
    }
    return result;
  }

  /**
   * Convert a reflection table from an array
   * @param array The array of reflections
   * @returns The reflection table
   */
  inline
  af::reflection_table reflection_table_from_array(af::const_ref<Reflection> array) {
    af::reflection_table result(array.size());
    for (std::size_t i = 0; i < array.size(); ++i) {
      detail::reflection_table_set_reflection(result, i, array[i]);
    }
    return result;
  }

}} // namespace dials::af

#endif // DIALS_ARRAY_FAMILY_REFLECTION_TABLE_H
