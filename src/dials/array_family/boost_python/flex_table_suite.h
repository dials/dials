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

#include <dxtbx/array_family/flex_table_suite.h>

namespace dials { namespace af { namespace boost_python { namespace flex_table_suite {

  using column_to_object_visitor [[deprecated(
    "column_to_object_visitor has moved to dxtbx/array_family/flex_table_suite.h")]] =
    dxtbx::af::flex_table_suite::column_to_object_visitor;
  using element_to_object_visitor [[deprecated(
    "element_to_object_visitor has moved to dxtbx/array_family/flex_table_suite.h")]] =
    dxtbx::af::flex_table_suite::element_to_object_visitor;
  using setitem_row_visitor [[deprecated(
    "setitem_row_visitor has moved to dxtbx/array_family/flex_table_suite.h")]] =
    dxtbx::af::flex_table_suite::setitem_row_visitor;
  using extend_column_visitor [[deprecated(
    "extend_column_visitor has moved to dxtbx/array_family/flex_table_suite.h")]] =
    dxtbx::af::flex_table_suite::extend_column_visitor;
  using update_column_visitor [[deprecated(
    "update_column_visitor has moved to dxtbx/array_family/flex_table_suite.h")]] =
    dxtbx::af::flex_table_suite::update_column_visitor;
  using copy_to_slice_visitor [[deprecated(
    "copy_to_slice_visitor has moved to dxtbx/array_family/flex_table_suite.h")]] =
    dxtbx::af::flex_table_suite::copy_to_slice_visitor;
  using copy_from_slice_visitor [[deprecated(
    "copy_from_slice_visitor has moved to dxtbx/array_family/flex_table_suite.h")]] =
    dxtbx::af::flex_table_suite::copy_from_slice_visitor;
  using copy_column_visitor [[deprecated(
    "copy_column_visitor has moved to dxtbx/array_family/flex_table_suite.h")]] =
    dxtbx::af::flex_table_suite::copy_column_visitor;
  using copy_from_indices_visitor [[deprecated(
    "copy_from_indices_visitor has moved to dxtbx/array_family/flex_table_suite.h")]] =
    dxtbx::af::flex_table_suite::copy_from_indices_visitor;
  using copy_to_indices_visitor [[deprecated(
    "copy_to_indices_visitor has moved to dxtbx/array_family/flex_table_suite.h")]] =
    dxtbx::af::flex_table_suite::copy_to_indices_visitor;
  using copy_to_indices_with_mask_visitor
    [[deprecated("copy_to_indices_with_mask_visitor has moved to "
                 "dxtbx/array_family/flex_table_suite.h")]] =
      dxtbx::af::flex_table_suite::copy_to_indices_with_mask_visitor;
  using reorder_visitor [[deprecated(
    "reorder_visitor has moved to dxtbx/array_family/flex_table_suite.h")]] =
    dxtbx::af::flex_table_suite::reorder_visitor;
  using compare_index
    [[deprecated("compare_index has moved to dxtbx/array_family/flex_table_suite.h")]] =
      dxtbx::af::flex_table_suite::compare_index;
  using sort_visitor
    [[deprecated("sort_visitor has moved to dxtbx/array_family/flex_table_suite.h")]] =
      dxtbx::af::flex_table_suite::sort_visitor;
  using remove_if_flag_visitor [[deprecated(
    "remove_if_flag_visitor has moved to dxtbx/array_family/flex_table_suite.h")]] =
    dxtbx::af::flex_table_suite::remove_if_flag_visitor;
  using make_flex_table [[deprecated(
    "make_flex_table has moved to dxtbx/array_family/flex_table_suite.h")]] =
    dxtbx::af::flex_table_suite::make_flex_table;
  using getitem_column [[deprecated(
    "getitem_column has moved to dxtbx/array_family/flex_table_suite.h")]] =
    dxtbx::af::flex_table_suite::getitem_column;
  using delitem_column [[deprecated(
    "delitem_column has moved to dxtbx/array_family/flex_table_suite.h")]] =
    dxtbx::af::flex_table_suite::delitem_column;
  using setitem_column [[deprecated(
    "setitem_column has moved to dxtbx/array_family/flex_table_suite.h")]] =
    dxtbx::af::flex_table_suite::setitem_column;
  using getitem_row
    [[deprecated("getitem_row has moved to dxtbx/array_family/flex_table_suite.h")]] =
      dxtbx::af::flex_table_suite::getitem_row;
  using delitem_row
    [[deprecated("delitem_row has moved to dxtbx/array_family/flex_table_suite.h")]] =
      dxtbx::af::flex_table_suite::delitem_row;
  using setitem_row
    [[deprecated("setitem_row has moved to dxtbx/array_family/flex_table_suite.h")]] =
      dxtbx::af::flex_table_suite::setitem_row;
  using getitem_slice
    [[deprecated("getitem_slice has moved to dxtbx/array_family/flex_table_suite.h")]] =
      dxtbx::af::flex_table_suite::getitem_slice;
  using remove_if_flag [[deprecated(
    "remove_if_flag has moved to dxtbx/array_family/flex_table_suite.h")]] =
    dxtbx::af::flex_table_suite::remove_if_flag;
  using delitem_slice
    [[deprecated("delitem_slice has moved to dxtbx/array_family/flex_table_suite.h")]] =
      dxtbx::af::flex_table_suite::delitem_slice;
  using setitem_slice
    [[deprecated("setitem_slice has moved to dxtbx/array_family/flex_table_suite.h")]] =
      dxtbx::af::flex_table_suite::setitem_slice;
  using has_key
    [[deprecated("has_key has moved to dxtbx/array_family/flex_table_suite.h")]] =
      dxtbx::af::flex_table_suite::has_key;
  using append
    [[deprecated("append has moved to dxtbx/array_family/flex_table_suite.h")]] =
      dxtbx::af::flex_table_suite::append;
  using insert
    [[deprecated("insert has moved to dxtbx/array_family/flex_table_suite.h")]] =
      dxtbx::af::flex_table_suite::insert;
  using update
    [[deprecated("update has moved to dxtbx/array_family/flex_table_suite.h")]] =
      dxtbx::af::flex_table_suite::update;
  using set_selected_rows_index [[deprecated(
    "set_selected_rows_index has moved to dxtbx/array_family/flex_table_suite.h")]] =
    dxtbx::af::flex_table_suite::set_selected_rows_index;
  using set_selected_rows_index_mask
    [[deprecated("set_selected_rows_index_mask has moved to "
                 "dxtbx/array_family/flex_table_suite.h")]] =
      dxtbx::af::flex_table_suite::set_selected_rows_index_mask;
  using set_selected_rows_flags [[deprecated(
    "set_selected_rows_flags has moved to dxtbx/array_family/flex_table_suite.h")]] =
    dxtbx::af::flex_table_suite::set_selected_rows_flags;
  using set_selected_cols_keys [[deprecated(
    "set_selected_cols_keys has moved to dxtbx/array_family/flex_table_suite.h")]] =
    dxtbx::af::flex_table_suite::set_selected_cols_keys;
  using set_selected_cols_tuple [[deprecated(
    "set_selected_cols_tuple has moved to dxtbx/array_family/flex_table_suite.h")]] =
    dxtbx::af::flex_table_suite::set_selected_cols_tuple;
  using del_selected_rows_flags [[deprecated(
    "del_selected_rows_flags has moved to dxtbx/array_family/flex_table_suite.h")]] =
    dxtbx::af::flex_table_suite::del_selected_rows_flags;
  using del_selected_rows_index [[deprecated(
    "del_selected_rows_index has moved to dxtbx/array_family/flex_table_suite.h")]] =
    dxtbx::af::flex_table_suite::del_selected_rows_index;
  using del_selected_cols_keys [[deprecated(
    "del_selected_cols_keys has moved to dxtbx/array_family/flex_table_suite.h")]] =
    dxtbx::af::flex_table_suite::del_selected_cols_keys;
  using del_selected_cols_tuple [[deprecated(
    "del_selected_cols_tuple has moved to dxtbx/array_family/flex_table_suite.h")]] =
    dxtbx::af::flex_table_suite::del_selected_cols_tuple;
  using type_appender
    [[deprecated("type_appender has moved to dxtbx/array_family/flex_table_suite.h")]] =
      dxtbx::af::flex_table_suite::type_appender;
  using types
    [[deprecated("types has moved to dxtbx/array_family/flex_table_suite.h")]] =
      dxtbx::af::flex_table_suite::types;
  using reorder
    [[deprecated("reorder has moved to dxtbx/array_family/flex_table_suite.h")]] =
      dxtbx::af::flex_table_suite::reorder;
  using sort [[deprecated("sort has moved to dxtbx/array_family/flex_table_suite.h")]] =
    dxtbx::af::flex_table_suite::sort;
  using copy [[deprecated("copy has moved to dxtbx/array_family/flex_table_suite.h")]] =
    dxtbx::af::flex_table_suite::copy;
  using deepcopy
    [[deprecated("deepcopy has moved to dxtbx/array_family/flex_table_suite.h")]] =
      dxtbx::af::flex_table_suite::deepcopy;
  using key_iterator
    [[deprecated("key_iterator has moved to dxtbx/array_family/flex_table_suite.h")]] =
      dxtbx::af::flex_table_suite::key_iterator;
  using column_iterator [[deprecated(
    "column_iterator has moved to dxtbx/array_family/flex_table_suite.h")]] =
    dxtbx::af::flex_table_suite::column_iterator;
  using row_iterator
    [[deprecated("row_iterator has moved to dxtbx/array_family/flex_table_suite.h")]] =
      dxtbx::af::flex_table_suite::row_iterator;
  using make_iterator
    [[deprecated("make_iterator has moved to dxtbx/array_family/flex_table_suite.h")]] =
      dxtbx::af::flex_table_suite::make_iterator;
  using flex_table_pickle_suite [[deprecated(
    "flex_table_pickle_suite has moved to dxtbx/array_family/flex_table_suite.h")]] =
    dxtbx::af::flex_table_suite::flex_table_pickle_suite;
  using setitem_column_generator [[deprecated(
    "setitem_column_generator has moved to dxtbx/array_family/flex_table_suite.h")]] =
    dxtbx::af::flex_table_suite::setitem_column_generator;
  using flex_table_wrapper [[deprecated(
    "flex_table_wrapper has moved to dxtbx/array_family/flex_table_suite.h")]] =
    dxtbx::af::flex_table_suite::flex_table_wrapper;

  // Reflection table specific
  using reflection_table_extend_identifiers
    [[deprecated("reflection_table_extend_identifers has moved to "
                 "dials/array_family/boost_python/reflection_table_suite.h")]] =
      dials::af::boost_python::reflection_table_suite::extend_identifiers;
  using extend [[deprecated(
    "extend has moved to dials/array_family/boost_python/reflection_table_suite.h")]] =
    dials::af::boost_python::reflection_table_suite::extend;
  using select_rows_index
    [[deprecated("select_rows_index has moved to "
                 "dials/array_family/boost_python/reflection_table_suite.h")]] =
      dials::af::boost_python::reflection_table_suite::select_rows_index;
  using select_rows_flags
    [[deprecated("select_rows_flags has moved to "
                 "dials/array_family/boost_python/reflection_table_suite.h")]] =
      dials::af::boost_python::reflection_table_suite::select_rows_flags;
  using select_cols_keys
    [[deprecated("select_cols_keys has moved to "
                 "dials/array_family/boost_python/reflection_table_suite.h")]] =
      dials::af::boost_python::reflection_table_suite::select_cols_keys;
  using select_cols_tuple
    [[deprecated("select_cols_tuple has moved to "
                 "dials/array_family/boost_python/reflection_table_suite.h")]] =
      dials::af::boost_python::reflection_table_suite::select_cols_tuple;
  using select_using_experiment
    [[deprecated("select_using_experiment has moved to "
                 "dials/array_family/boost_python/reflection_table_suite.h")]] =
      dials::af::boost_python::reflection_table_suite::select_using_experiment;
  using select_using_experiments
    [[deprecated("select_using_experiments has moved to "
                 "dials/array_family/boost_python/reflection_table_suite.h")]] =
      dials::af::boost_python::reflection_table_suite::select_using_experiments;

}}}}  // namespace dials::af::boost_python::flex_table_suite

#endif  // DIALS_FRAMEWORK_TABLE_BOOST_PYTHON_FLEX_TABLE_SUITE_H
