/*
 * flex_table_suite.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */

#ifndef DIALS_FRAMEWORK_TABLE_BOOST_PYTHON_FLEX_TABLE_SUITE_H
#define DIALS_FRAMEWORK_TABLE_BOOST_PYTHON_FLEX_TABLE_SUITE_H

#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <boost/python/slice.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <dxtbx/array_family/flex_table_suite.h>
#include <dials/array_family/boost_python/reflection_table_suite.h>

namespace dials { namespace af { namespace boost_python { namespace flex_table_suite {

  using namespace boost::python;

  using column_to_object_visitor [[deprecated(
    "column_to_object_visitor has moved to dxtbx/array_family/flex_table_suite.h")]] =
    dxtbx::af::flex_table_suite::column_to_object_visitor;

  using element_to_object_visitor [[deprecated(
    "element_to_object_visitor has moved to dxtbx/array_family/flex_table_suite.h")]] =
    dxtbx::af::flex_table_suite::element_to_object_visitor;

  using setitem_row_visitor [[deprecated(
    "setitem_row_visitor has moved to dxtbx/array_family/flex_table_suite.h")]] =
    dxtbx::af::flex_table_suite::setitem_row_visitor;

  template <typename T>
  using extend_column_visitor [[deprecated(
    "extend_column_visitor has moved to dxtbx/array_family/flex_table_suite.h")]] =
    dxtbx::af::flex_table_suite::extend_column_visitor<T>;

  template <typename T>
  using update_column_visitor [[deprecated(
    "update_column_visitor has moved to dxtbx/array_family/flex_table_suite.h")]] =
    dxtbx::af::flex_table_suite::update_column_visitor<T>;

  template <typename T>
  using copy_to_slice_visitor [[deprecated(
    "copy_to_slice_visitor has moved to dxtbx/array_family/flex_table_suite.h")]] =
    dxtbx::af::flex_table_suite::copy_to_slice_visitor<T>;

  template <typename T>
  using copy_from_slice_visitor [[deprecated(
    "copy_from_slice_visitor has moved to dxtbx/array_family/flex_table_suite.h")]] =
    dxtbx::af::flex_table_suite::copy_from_slice_visitor<T>;

  template <typename T>
  using copy_column_visitor [[deprecated(
    "copy_column_visitor has moved to dxtbx/array_family/flex_table_suite.h")]] =
    dxtbx::af::flex_table_suite::copy_column_visitor<T>;

  template <typename T>
  using copy_from_indices_visitor [[deprecated(
    "copy_from_indices_visitor has moved to dxtbx/array_family/flex_table_suite.h")]] =
    dxtbx::af::flex_table_suite::copy_from_indices_visitor<T>;

  template <typename T>
  using copy_to_indices_visitor [[deprecated(
    "copy_to_indices_visitor has moved to dxtbx/array_family/flex_table_suite.h")]] =
    dxtbx::af::flex_table_suite::copy_to_indices_visitor<T>;

  template <typename T>
  using copy_to_indices_with_mask_visitor
    [[deprecated("copy_to_indices_with_mask_visitor has moved to "
                 "dxtbx/array_family/flex_table_suite.h")]] =
      dxtbx::af::flex_table_suite::copy_to_indices_with_mask_visitor<T>;

  using reorder_visitor [[deprecated(
    "reorder_visitor has moved to dxtbx/array_family/flex_table_suite.h")]] =
    dxtbx::af::flex_table_suite::reorder_visitor;

  template <typename T>
  using compare_index
    [[deprecated("compare_index has moved to dxtbx/array_family/flex_table_suite.h")]] =
      dxtbx::af::flex_table_suite::compare_index<T>;

  using sort_visitor
    [[deprecated("sort_visitor has moved to dxtbx/array_family/flex_table_suite.h")]] =
      dxtbx::af::flex_table_suite::sort_visitor;

  using remove_if_flag_visitor [[deprecated(
    "remove_if_flag_visitor has moved to dxtbx/array_family/flex_table_suite.h")]] =
    dxtbx::af::flex_table_suite::remove_if_flag_visitor;

  template <typename T>
  [[deprecated(
    "make_flex_table has moved to dxtbx/array_family/flex_table_suite.h")]] T *
  make_flex_table(list columns) {
    return dxtbx::af::flex_table_suite::make_flex_table<T>(columns);
  }

  template <typename T>
  [[deprecated(
    "getitem_column has moved to dxtbx/array_family/flex_table_suite.h")]] object
  getitem_column(T &self, const typename T::key_type &key) {
    return dxtbx::af::flex_table_suite::getitem_column(self, key);
  }

  template <typename T>
  [[deprecated(
    "delitem_column has moved to dxtbx/array_family/flex_table_suite.h")]] void
  delitem_column(T &self, const typename T::key_type &key) {
    dxtbx::af::flex_table_suite::delitem_column(self, key);
  }

  template <typename T, typename U>
  [[deprecated(
    "setitem_column has moved to dxtbx/array_family/flex_table_suite.h")]] void
  setitem_column(T &self,
                 const typename T::key_type &key,
                 const scitbx::af::const_ref<U> &data) {
    dxtbx::af::flex_table_suite::setitem_column(self, key, data);
  }

  template <typename T>
  [[deprecated("getitem_row has moved to dxtbx/array_family/flex_table_suite.h")]] dict
  getitem_row(const T &self, typename T::size_type n) {
    dxtbx::af::flex_table_suite::getitem_row(self, n);
  }

  template <typename T>
  [[deprecated("delitem_row has moved to dxtbx/array_family/flex_table_suite.h")]] void
  delitem_row(T &self, typename T::size_type n) {
    dxtbx::af::flex_table_suite::delitem_row(self, n);
  }

  template <typename T>
  [[deprecated("setitem_row has moved to dxtbx/array_family/flex_table_suite.h")]] void
  setitem_row(T &self, typename T::size_type n, dict row) {
    dxtbx::af::flex_table_suite::setitem_row(self, n, row);
  }

  template <typename T>
  [[deprecated(
    "remove_if_flag has moved to dxtbx/array_family/flex_table_suite.h")]] void
  remove_if_flag(T &self, const scitbx::af::const_ref<bool> &flags) {
    dxtbx::af::flex_table_suite::remove_if_flag(self, flags);
  }

  template <typename T>
  [[deprecated(
    "delitem_slice has moved to dxtbx/array_family/flex_table_suite.h")]] void
  delitem_slice(T &self, slice s) {
    dxtbx::af::flex_table_suite::delitem_slice(self, s);
  }

  template <typename T>
  [[deprecated(
    "setitem_slice has moved to dxtbx/array_family/flex_table_suite.h")]] void
  setitem_slice(T &self, slice s, const T &other) {
    dxtbx::af::flex_table_suite::setitem_slice(self, s, other);
  }

  template <typename T>
  [[deprecated("has_key has moved to dxtbx/array_family/flex_table_suite.h")]] bool
  has_key(const T &self, const typename T::key_type &key) {
    return dxtbx::af::flex_table_suite::has_key(self, key);
  }

  template <typename T>
  [[deprecated("append has moved to dxtbx/array_family/flex_table_suite.h")]] void
  append(T &self, dict row) {
    dxtbx::af::flex_table_suite::append(self, row);
  }

  template <typename T>
  [[deprecated("insert has moved to dxtbx/array_family/flex_table_suite.h")]] void
  insert(T &self, typename T::size_type n, dict row) {
    dxtbx::af::flex_table_suite::insert(self, n, row);
  }

  template <typename T>
  [[deprecated("update has moved to dxtbx/array_family/flex_table_suite.h")]] void
  update(T &self, const T &other) {
    dxtbx::af::flex_table_suite::update(self, other);
  }

  template <typename T>
  [[deprecated(
    "set_selected_rows_index has moved to dxtbx/array_family/flex_table_suite.h")]] void
  set_selected_rows_index(T &self,
                          const scitbx::af::const_ref<std::size_t> &index,
                          const T &other) {
    dxtbx::af::flex_table_suite::set_selected_rows_index(self, index, other);
  }

  template <typename T>
  [[deprecated(
    "set_selected_rows_index_mask has moved to "
    "dxtbx/array_family/flex_table_suite.h")]] void
  set_selected_rows_index_mask(T &self,
                               const scitbx::af::const_ref<std::size_t> &index,
                               const scitbx::af::const_ref<bool> &mask,
                               const T &other) {
    dxtbx::af::flex_table_suite::set_selected_rows_index_mask(self, index, mask, other);
  }

  template <typename T>
  [[deprecated(
    "set_selected_rows_flags has moved to dxtbx/array_family/flex_table_suite.h")]] void
  set_selected_rows_flags(T &self,
                          const scitbx::af::const_ref<bool> &flags,
                          const T &other) {
    dxtbx::af::flex_table_suite::set_selected_rows_flags(self, flags, other);
  }

  template <typename T>
  [[deprecated(
    "set_selected_cols_keys has moved to dxtbx/array_family/flex_table_suite.h")]] void
  set_selected_cols_keys(T &self,
                         const scitbx::af::const_ref<std::string> &keys,
                         const T &other) {
    dxtbx::af::flex_table_suite::set_selected_cols_keys(self, keys, other);
  }

  template <typename T>
  [[deprecated(
    "set_selected_cols_tuple has moved to dxtbx/array_family/flex_table_suite.h")]] void
  set_selected_cols_tuple(T &self, boost::python::tuple keys, const T &other) {
    dxtbx::af::flex_table_suite::set_selected_cols_tuple(keys, other);
  }

  template <typename T>
  [[deprecated(
    "del_selected_rows_flags has moved to dxtbx/array_family/flex_table_suite.h")]] void
  del_selected_rows_flags(T &self, const scitbx::af::const_ref<bool> &flags) {
    dxtbx::af::flex_table_suite::del_selected_rows_flags(self, flags);
  }

  template <typename T>
  [[deprecated(
    "del_selected_rows_index has moved to dxtbx/array_family/flex_table_suite.h")]] void
  del_selected_rows_index(T &self, const scitbx::af::const_ref<std::size_t> &index) {
    dxtbx::af::flex_table_suite::del_selected_rows_index(self, index);
  }

  template <typename T>
  [[deprecated(
    "del_selected_cols_keys has moved to dxtbx/array_family/flex_table_suite.h")]] void
  del_selected_cols_keys(T &self, const scitbx::af::const_ref<std::string> &keys) {
    dxtbx::af::flex_table_suite::del_selected_cols_keys(self, keys);
  }

  template <typename T>
  [[deprecated(
    "del_selected_cols_tuple has moved to dxtbx/array_family/flex_table_suite.h")]] void
  del_selected_cols_tuple(T &self, boost::python::tuple keys) {
    dxtbx::af::flex_table_suite::del_selected_cols_tuple(self, keys);
  }

  using type_appender
    [[deprecated("type_appender has moved to dxtbx/array_family/flex_table_suite.h")]] =
      dxtbx::af::flex_table_suite::type_appender;

  template <typename T>
  [[deprecated("types has moved to dxtbx/array_family/flex_table_suite.h")]] list types(
    const T &self) {
    return dxtbx::af::flex_table_suite::types(self);
  }

  template <typename T>
  [[deprecated("reorder has moved to dxtbx/array_family/flex_table_suite.h")]] void
  reorder(T &self, const scitbx::af::const_ref<std::size_t> &index) {
    dxtbx::af::flex_table_suite::reorder(self, index);
  }

  template <typename T>
  [[deprecated("sort has moved to dxtbx/array_family/flex_table_suite.h")]] void
  sort(T &self, typename T::key_type key, bool reverse) {
    dxtbx::af::flex_table_suite::sort(self, key, reverse);
  }

  template <typename T>
  [[deprecated("copy has moved to dxtbx/array_family/flex_table_suite.h")]] T copy(
    const T &self) {
    return dxtbx::af::flex_table_suite::copy(self);
  }

  template <typename T>
  [[deprecated("deepcopy has moved to dxtbx/array_family/flex_table_suite.h")]] T
  deepcopy(const T &self, dict obj) {
    return dxtbx::af::flex_table_suite::deepcopy(self, obj);
  }

  template <typename T>
  using key_iterator
    [[deprecated("key_iterator has moved to dxtbx/array_family/flex_table_suite.h")]] =
      dxtbx::af::flex_table_suite::key_iterator<T>;

  template <typename T>
  using column_iterator [[deprecated(
    "column_iterator has moved to dxtbx/array_family/flex_table_suite.h")]] =
    dxtbx::af::flex_table_suite::column_iterator<T>;

  template <typename T>
  using row_iterator
    [[deprecated("row_iterator has moved to dxtbx/array_family/flex_table_suite.h")]] =
      dxtbx::af::flex_table_suite::row_iterator<T>;

  template <typename T>
  using make_iterator
    [[deprecated("make_iterator has moved to dxtbx/array_family/flex_table_suite.h")]] =
      dxtbx::af::flex_table_suite::make_iterator<T>;

  template <typename T>
  using flex_table_pickle_suite [[deprecated(
    "flex_table_pickle_suite has moved to dxtbx/array_family/flex_table_suite.h")]] =
    dxtbx::af::flex_table_suite::flex_table_pickle_suite<T>;

  template <typename T>
  using setitem_column_generator [[deprecated(
    "setitem_column_generator has moved to dxtbx/array_family/flex_table_suite.h")]] =
    dxtbx::af::flex_table_suite::setitem_column_generator<T>;

  template <typename T>
  using flex_table_wrapper [[deprecated(
    "flex_table_wrapper has moved to dxtbx/array_family/flex_table_suite.h")]] =
    dxtbx::af::flex_table_suite::flex_table_wrapper<T>;

  // Reflection table specific
  template <typename T>
  [[deprecated(
    "reflection_table_extend_identifers has moved to "
    "dials/array_family/boost_python/reflection_table_suite.h")]] void
  reflection_table_extend_identifiers(T &self, const T &other) {
    dials::af::boost_python::reflection_table_suite::extend_identifiers(self, other);
  }

  template <typename T>
  [[deprecated(
    "extend has moved to "
    "dials/array_family/boost_python/reflection_table_suite.h")]] void
  extend(T &self, const T &other) {
    return dials::af::boost_python::reflection_table_suite::extend(self, other);
  }

  template <typename T>
  [[deprecated(
    "select_rows_index has moved to "
    "dials/array_family/boost_python/reflection_table_suite.h")]] T
  select_rows_index(const T &self, const scitbx::af::const_ref<std::size_t> &index) {
    return dials::af::boost_python::reflection_table_suite::select_rows_index(self,
                                                                              index);
  }

  template <typename T>
  [[deprecated(
    "select_rows_flags has moved to "
    "dials/array_family/boost_python/reflection_table_suite.h")]] T
  select_rows_flags(const T &self, const scitbx::af::const_ref<bool> &flags) {
    return dials::af::boost_python::reflection_table_suite::select_rows_flags(self,
                                                                              flags);
  }

  template <typename T>
  [[deprecated(
    "select_cols_keys has moved to "
    "dials/array_family/boost_python/reflection_table_suite.h")]] T
  select_cols_keys(const T &self, const scitbx::af::const_ref<std::string> &keys) {
    return dials::af::boost_python::reflection_table_suite::select_cols_keys(self,
                                                                             keys);
  }

  template <typename T>
  [[deprecated(
    "select_cols_tuple has moved to "
    "dials/array_family/boost_python/reflection_table_suite.h")]] T
  select_cols_tuple(const T &self, boost::python::tuple keys) {
    return dials::af::boost_python::reflection_table_suite::select_cols_tuple(self,
                                                                              keys);
  }

  template <typename T>
  [[deprecated(
    "select_using_experiment has moved to "
    "dials/array_family/boost_python/reflection_table_suite.h")]] T
  select_using_experiment(T &self, dxtbx::model::Experiment expt) {
    return dials::af::boost_python::reflection_table_suite::select_using_experiment(
      self, expt);
  }

  template <typename T>
  [[deprecated(
    "select_using_experiments has moved to "
    "dials/array_family/boost_python/reflection_table_suite.h")]] T
  select_using_experiments(T &self, dxtbx::model::ExperimentList expts) {
    return dials::af::boost_python::reflection_table_suite::select_using_experiments(
      self, expts);
  }

  template <typename T>
  [[deprecated(
    "getitem_slice has moved to "
    "dials/array_family/boost_python/reflection_table_suite.h")]] T
  getitem_slice(const T &self, slice s) {
    return dials::af::boost_python::reflection_table_suite::getitem_slice(self, s);
  }

}}}}  // namespace dials::af::boost_python::flex_table_suite

#endif  // DIALS_FRAMEWORK_TABLE_BOOST_PYTHON_FLEX_TABLE_SUITE_H
