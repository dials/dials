/*
 * column_table.cc
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <boost/python/slice.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/mpl/for_each.hpp>
#include <string>
#include <iterator>
#include <iostream>
#include <sstream>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/boost_python/slice.h>
#include <dials/framework/table/column_table.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/error.h>

namespace dials { namespace framework { namespace boost_python {
namespace column_table_suite {

  using namespace boost::python;

  struct column_data_wrapper {
    std::string name_;

    column_data_wrapper(std::string name)
      : name_(name) {}

    template <typename T>
    void operator()(const T &x) {
      std::string name = name_ + typeid(x).name();
      typedef T column_data_type;
      class_<column_data_type>(name.c_str())
        .def(vector_indexing_suite<column_data_type>())
        .def("resize", &column_data_type::resize);
    }
  };



  template <typename T, typename U>
  void column_table_set_data(
      T &table,
      const typename T::key_type &key,
      const af::const_ref<U> &a) {
    table.erase(key);
    column_data<U> col = table[key];
    col.resize(a.size());
    std::copy(a.begin(), a.end(), col.begin());
  }

  template <typename T>
  struct column_table_set_data_wrapper {
    class_<T> table_class_;
    column_table_set_data_wrapper(class_<T> table_class)
      : table_class_(table_class) {}
    template <typename U>
    void operator()(const U &x) {
      table_class_.def("__setitem__",
        &column_table_set_data<T, typename U::value_type>,
          return_internal_reference<>());
    }
  };

  struct column_data_to_object :
      public boost::static_visitor<boost::python::object> {
    template <typename T>
    boost::python::object operator () (T &col) {
      return boost::python::object(col);
    }
  };

  template <typename ColumnTable>
  boost::python::object column_table_get_data(
      ColumnTable &table,
      const typename ColumnTable::key_type &key) {
    typedef typename ColumnTable::mapped_type mapped_type;
    mapped_type col = table[key];
    DIALS_ASSERT(!col.empty());
    column_data_to_object to_object;
    return col.apply_visitor(to_object);
  }

  template <typename ColumnTable>
  bool column_table_has_key(
      const ColumnTable &table,
      const typename ColumnTable::key_type &key) {
    return table.count(key) == 1;
  }

  template <typename ColumnTable>
  boost::python::list column_table_items(const ColumnTable &table) {
    boost::python::list result;
    column_data_to_object to_object;
    for (typename ColumnTable::const_iterator it = table.begin();
        it != table.end(); ++it) {
      result.append(boost::python::make_tuple(
        it->first,
        it->second.apply_visitor(to_object)));
    }
    return result;
  }

  template <typename ColumnTable>
  boost::python::list column_table_keys(const ColumnTable &table) {
    boost::python::list result;
    for (typename ColumnTable::const_iterator it = table.begin();
        it != table.end(); ++it) {
      result.append(it->first);
    }
    return result;
  }

  template <typename ColumnTable>
  boost::python::list column_table_values(const ColumnTable &table) {
    boost::python::list result;
    column_data_to_object to_object;
    for (typename ColumnTable::const_iterator it = table.begin();
        it != table.end(); ++it) {
      result.append(it->second.apply_visitor(to_object));
    }
    return result;
  }

  template <typename ColumnTable>
  void column_table_update(ColumnTable &self, const ColumnTable &other) {

  }

  template <typename ColumnTable>
  class column_table_iterkeys_proxy {
  public:

    typedef typename ColumnTable::const_iterator base_iterator;

    typedef ptrdiff_t difference_type;
    typedef std::forward_iterator_tag iterator_category;
    typedef typename ColumnTable::key_type value_type;
    typedef const value_type *pointer;
    typedef const value_type &reference;

    column_table_iterkeys_proxy(base_iterator it)
      : it_(it) {}

    reference operator*() {
      return it_->first;
    }

    column_table_iterkeys_proxy& operator++() {
      ++it_;
      return *this;
    }

    column_table_iterkeys_proxy operator++(int) {
      column_table_iterkeys_proxy result(*this);
      ++(*this);
      return result;
    }

    bool operator==(const column_table_iterkeys_proxy& rhs) const {
      return it_ == rhs.it_;
    }

    bool operator!=(const column_table_iterkeys_proxy& rhs) const {
      return !(*this == rhs);
    }

  private:
    base_iterator it_;
  };

  template <typename ColumnTable>
  column_table_iterkeys_proxy<ColumnTable> column_table_iterkeys_begin(
      const ColumnTable &table) {
    return column_table_iterkeys_proxy<ColumnTable>(table.begin());
  }

  template <typename ColumnTable>
  column_table_iterkeys_proxy<ColumnTable> column_table_iterkeys_end(
      const ColumnTable &table) {
    return column_table_iterkeys_proxy<ColumnTable>(table.end());
  }

  template <typename ColumnTable>
  class column_table_itervalues_proxy {
  public:

    typedef typename ColumnTable::const_iterator base_iterator;

    typedef ptrdiff_t difference_type;
    typedef std::forward_iterator_tag iterator_category;
    typedef boost::python::object value_type;
    typedef const value_type *pointer;
    typedef const value_type reference;

    column_table_itervalues_proxy(base_iterator it)
      : it_(it) {}

    reference operator*() {
      column_data_to_object to_object;
      return it_->second.apply_visitor(to_object);
    }

    column_table_itervalues_proxy& operator++() {
      ++it_;
      return *this;
    }

    column_table_itervalues_proxy operator++(int) {
      column_table_itervalues_proxy result(*this);
      ++(*this);
      return result;
    }

    bool operator==(const column_table_itervalues_proxy& rhs) const {
      return it_ == rhs.it_;
    }

    bool operator!=(const column_table_itervalues_proxy& rhs) const {
      return !(*this == rhs);
    }

  private:
    base_iterator it_;
  };


  template <typename ColumnTable>
  column_table_itervalues_proxy<ColumnTable> column_table_itervalues_begin(
      const ColumnTable &table) {
    return column_table_itervalues_proxy<ColumnTable>(table.begin());
  }

  template <typename ColumnTable>
  column_table_itervalues_proxy<ColumnTable> column_table_itervalues_end(
      const ColumnTable &table) {
    return column_table_itervalues_proxy<ColumnTable>(table.end());
  }


  template <typename ColumnTable>
  class column_table_iteritems_proxy {
  public:

    typedef typename ColumnTable::const_iterator base_iterator;

    typedef ptrdiff_t difference_type;
    typedef std::forward_iterator_tag iterator_category;
    typedef boost::python::object value_type;
    typedef const value_type *pointer;
    typedef const value_type reference;

    column_table_iteritems_proxy(base_iterator it)
      : it_(it) {}

    reference operator*() {
      column_data_to_object to_object;
      return boost::python::make_tuple(
        it_->first,
        it_->second.apply_visitor(to_object));
    }

    column_table_iteritems_proxy& operator++() {
      ++it_;
      return *this;
    }

    column_table_iteritems_proxy operator++(int) {
      column_table_iteritems_proxy result(*this);
      ++(*this);
      return result;
    }

    bool operator==(const column_table_iteritems_proxy& rhs) const {
      return it_ == rhs.it_;
    }

    bool operator!=(const column_table_iteritems_proxy& rhs) const {
      return !(*this == rhs);
    }

  private:
    base_iterator it_;
  };

  template <typename ColumnTable>
  class column_table_iterrows_proxy {
  public:

    typedef ptrdiff_t difference_type;
    typedef std::forward_iterator_tag iterator_category;
    typedef boost::python::object value_type;
    typedef const value_type *pointer;
    typedef const value_type reference;

    column_table_iterrows_proxy(
        typename ColumnTable::const_iterator first,
        typename ColumnTable::const_iterator last,
        std::size_t index)
      : index_(index) {
      for (; first != last; ++first) {
        header_.push_back(first->first);
        columns_.push_back(first->second);
      }
    }

    reference operator*() {
      boost::python::dict result;
      column_data_to_object to_object;
      for (std::size_t i = 0; i < columns_.size(); ++i) {
        result[header_[i]] = columns_[i].apply_visitor(to_object)[index_];
      }
      return result;
    }

    column_table_iterrows_proxy& operator++() {
      ++index_;
      return *this;
    }

    column_table_iterrows_proxy operator++(int) {
      column_table_iterrows_proxy result(*this);
      ++(*this);
      return result;
    }

    bool operator==(const column_table_iterrows_proxy& rhs) const {
      return index_ == rhs.index_;
    }

    bool operator!=(const column_table_iterrows_proxy& rhs) const {
      return !(*this == rhs);
    }

  private:
    std::size_t index_;
    std::vector<typename ColumnTable::mapped_type> columns_;
    std::vector<std::string> header_;
  };

  template <typename ColumnTable>
  column_table_iteritems_proxy<ColumnTable> column_table_iteritems_begin(
      const ColumnTable &table) {
    return column_table_iteritems_proxy<ColumnTable>(table.begin());
  }

  template <typename ColumnTable>
  column_table_iteritems_proxy<ColumnTable> column_table_iteritems_end(
      const ColumnTable &table) {
    return column_table_iteritems_proxy<ColumnTable>(table.end());
  }


  template <typename ColumnTable>
  column_table_iterrows_proxy<ColumnTable> column_table_iterrows_begin(
      const ColumnTable &table) {
    return column_table_iterrows_proxy<ColumnTable>(table.begin(), table.end(), 0);
  }

  template <typename ColumnTable>
  column_table_iterrows_proxy<ColumnTable> column_table_iterrows_end(
      const ColumnTable &table) {
    return column_table_iterrows_proxy<ColumnTable>(table.begin(), table.end(), table.nrows());
  }

  template <typename ColumnTable>
  boost::python::object column_table_get_row_data(
      ColumnTable &table, std::size_t index) {
    boost::python::dict result;
    column_data_to_object to_object;
    for (typename ColumnTable::iterator it = table.begin(); it != table.end(); ++it) {
      result[it->first] = it->second.apply_visitor(to_object)[index];
    }
    return result;
  }

  template <typename ColumnTable>
  struct copy_column_data_slice :
      public boost::static_visitor<void> {

    ColumnTable &table;
    typename ColumnTable::key_type name;
    scitbx::boost_python::adapted_slice slice;

    copy_column_data_slice(
          ColumnTable &table_,
          typename ColumnTable::key_type name_,
          scitbx::boost_python::adapted_slice slice_)
      : table(table_),
        name(name_),
        slice(slice_) {}

    template <typename T>
    void operator () (const T &col) {
      T c = table[name];
      for (std::size_t i = 0, j = slice.start; i < table.nrows(); ++i, j += slice.step) {
        c[i] = col[j];
      }
    }
  };

  template <typename ColumnTable>
  struct copy_column_data_from_slice :
      public boost::static_visitor<void> {

    ColumnTable &table;
    typename ColumnTable::key_type name;
    scitbx::boost_python::adapted_slice slice;
    std::size_t num;

    copy_column_data_from_slice(
          ColumnTable &table_,
          typename ColumnTable::key_type name_,
          scitbx::boost_python::adapted_slice slice_,
          std::size_t num_)
      : table(table_),
        name(name_),
        slice(slice_),
        num(num_) {}

    template <typename T>
    void operator () (const T &col) {
      T c = table[name];
      for (std::size_t i = 0, j = slice.start; i < num; ++i, j += slice.step) {
        c[j] = col[i];
      }
    }
  };

  template <typename ColumnTable>
  ColumnTable column_table_get_slice(const ColumnTable &table,
      boost::python::slice slice) {

    scitbx::boost_python::adapted_slice aslice(slice, table.nrows());

    ColumnTable result(aslice.size);
    for (typename ColumnTable::const_iterator it = table.begin(); it != table.end(); ++it) {
      copy_column_data_slice <ColumnTable> copy_slice(result, it->first, aslice);
      it->second.apply_visitor(copy_slice);
    }
    return result;
  }

  template <typename ColumnTable>
  void column_table_set_slice(ColumnTable &table,
      boost::python::slice slice, const ColumnTable &other) {

    scitbx::boost_python::adapted_slice aslice(slice, table.nrows());
    for (typename ColumnTable::const_iterator it = other.begin(); it != other.end(); ++it) {
      copy_column_data_from_slice <ColumnTable> copy_slice(table, it->first, aslice, other.nrows());
      it->second.apply_visitor(copy_slice);
    }
  }


  struct set_row_data :
      public boost::static_visitor<void> {
    std::size_t index;
    boost::python::object item;
    set_row_data(std::size_t index_, boost::python::object item_)
      : index(index_),
        item(item_) {}

    template <typename T>
    void operator () (T &col) {
      col[index] = boost::python::extract<typename T::value_type>(item);
    }
  };

  template <typename ColumnTable>
  void column_table_set_row(ColumnTable &table,
      std::size_t index, boost::python::dict row) {
    for (typename ColumnTable::iterator it = table.begin(); it != table.end(); ++it) {
      if (row.has_key(it->first)) {
        set_row_data set_item(index, row[it->first]);
        it->second.apply_visitor(set_item);
      }
    }
  }

  template <typename ColumnTable>
  ColumnTable* make_column_table(boost::python::object columns) {
    ColumnTable *table = new ColumnTable();
    boost::python::object table_obj(table);
    for (std::size_t i = 0; i < boost::python::len(columns); ++i) {
      table_obj[columns[i][0]] = columns[i][1];
    }
    return table;
  }

  template <typename ColumnTable>
  void append(ColumnTable &table, boost::python::dict row) {
    table.resize(table.nrows() + 1);
    column_table_set_row(table, table.nrows() - 1, row);
  }

  template <typename ColumnTable>
  void insert(ColumnTable &table, std::size_t index, boost::python::dict row) {
    table.insert(index);
    column_table_set_row(table, index, row);
  }

  template <typename ColumnTable>
  struct copy_column_data :
      public boost::static_visitor<void> {

    ColumnTable &table;
    typename ColumnTable::key_type name;
    std::size_t na, nb;

    copy_column_data(
          ColumnTable &table_,
          typename ColumnTable::key_type name_,
          std::size_t na_, std::size_t nb_)
      : table(table_),
        name(name_),
        na(na_),
        nb(nb_) {}

    template <typename T>
    void operator () (const T &col) {
      T c = table[name];
      for (std::size_t i = 0; i < nb; ++i) {
        c[na + i] = col[i];
      }
    }
  };

  template <typename ColumnTable>
  void extend(ColumnTable &table, const ColumnTable &other) {
    std::size_t na = table.nrows();
    std::size_t nb = other.nrows();
    table.resize(na + nb);
    for (typename ColumnTable::const_iterator it = other.begin(); it != other.end(); ++it) {
      copy_column_data <ColumnTable> copy_data(table, it->first, na, nb);
      it->second.apply_visitor(copy_data);
    }
  }

  template <typename ColumnTable>
  struct add_column_data :
      public boost::static_visitor<void> {

    ColumnTable &table;
    typename ColumnTable::key_type name;

    add_column_data(
          ColumnTable &table_,
          typename ColumnTable::key_type name_)
      : table(table_),
        name(name_) {}

    template <typename T>
    void operator () (const T &col) {
      if (table.count(name)) {
        table.erase(name);
      }
      T c = table[name];
      for (std::size_t i = 0; i < col.size(); ++i) {
        c[i] = col[i];
      }
    }
  };

  template <typename ColumnTable>
  void update(ColumnTable &table, const ColumnTable &other) {
    for (typename ColumnTable::const_iterator it = other.begin(); it != other.end(); ++it) {
      add_column_data<ColumnTable> add_column(table, it->first);
      it->second.apply_visitor(add_column);
    }
  }

  struct reorder_column_data :
      public boost::static_visitor<void> {

    af::const_ref<std::size_t> index;

    reorder_column_data(const af::const_ref<std::size_t> &index_)
      : index(index_) {}

    template <typename T>
    void operator () (T col) {
      std::vector<typename T::value_type> temp(col.begin(), col.end());
      for (std::size_t i = 0; i < index.size(); ++i) {
        col[i] = temp[index[i]];
      }
    }
  };

  template <typename ColumnTable>
  void reorder_table(ColumnTable &table, const af::const_ref<std::size_t> &index) {
    reorder_column_data rec(index);
    for (typename ColumnTable::const_iterator it = table.begin(); it != table.end(); ++it) {
      it->second.apply_visitor(rec);
    }
  }

  template <typename T>
  struct compare_index {
    const T &v_;
    compare_index(const T &v)
      : v_(v) {}
    bool operator()(std::size_t a, std::size_t b) {
      return v_[a] < v_[b];
    }
  };

  struct sort_column_data :
      public boost::static_visitor<void> {
    af::ref<std::size_t> index;
    sort_column_data(af::ref<std::size_t> index_)
      : index(index_) {}

    template <typename T>
    void operator () (const T col) {
      for (std::size_t i = 0; i < index.size(); ++i) {
        index[i] = i;
      }
      std::sort(index.begin(), index.end(), compare_index<T>(col));
    }
  };

  template <typename ColumnTable>
  void sort_table(ColumnTable &table, std::string column, bool reverse) {
    af::shared<std::size_t> index(table.nrows());
    sort_column_data sc(index.ref());
    typename ColumnTable::mapped_type col = table[column];
    col.apply_visitor(sc);
    if (reverse) {
      std::reverse(index.begin(), index.end());
    }
    reorder_table(table, index.const_ref());
  }

  struct type_printer {
    boost::python::list type_list;
    type_printer(boost::python::list type_list_)
      : type_list(type_list_) {}

    template <typename U>
    void operator()(U x) {
      typedef typename U::value_type type;
      type a;
      boost::python::object obj(a);
      boost::python::object tt(boost::python::handle<>(PyObject_Type(obj.ptr())));
      type_list.append(tt);
    }
  };

  template <typename ColumnTable>
  boost::python::list types(const ColumnTable &table) {
    boost::python::list result;
    typedef typename ColumnTable::mapped_type::types column_types;
    boost::mpl::for_each<column_types>(type_printer(result));
    return result;
  }

  template <typename column_types>
  void column_table_wrapper(const char *name) {

    boost::mpl::for_each<typename column_types::types>(
      column_data_wrapper(name));

    typedef column_table<column_types> column_table_type;

    class_<column_table_type> column_table_class(name);
    column_table_class
      .def(init<std::size_t>())
      .def("__init__", make_constructor(&make_column_table<column_table_type>))
      .def("clear", &column_table_type::clear)
      .def("has_key", &column_table_has_key<column_table_type>)
      .def("items", &column_table_items<column_table_type>)
      .def("keys", &column_table_keys<column_table_type>)
      .def("values", &column_table_values<column_table_type>)
      .def("iteritems", range(
        &column_table_iteritems_begin<column_table_type>,
        &column_table_iteritems_end<column_table_type>))
      .def("iterkeys", range(
        &column_table_iterkeys_begin<column_table_type>,
        &column_table_iterkeys_end<column_table_type>))
      .def("itervalues", range(
        &column_table_itervalues_begin<column_table_type>,
        &column_table_itervalues_end<column_table_type>))
      .def("iterrows", range(
        &column_table_iterrows_begin<column_table_type>,
        &column_table_iterrows_end<column_table_type>))
      .def("update", &column_table_update<column_table_type>)
      //.def("erase", &column_table_type::erase)
      .def("empty", &column_table_type::empty)
      .def("nrows", &column_table_type::nrows)
      .def("ncols", &column_table_type::ncols)
      .def("__len__", &column_table_type::size)
      .def("__contains__", &column_table_has_key<column_table_type>)
      .def("__getitem__", &column_table_get_data<column_table_type>)
      .def("__getitem__", &column_table_get_row_data<column_table_type>)
      .def("__getitem__", &column_table_get_slice<column_table_type>)
      .def("__setitem__", &column_table_set_slice<column_table_type>)
      .def("__setitem__", &column_table_set_row<column_table_type>)
      .def("append", &append<column_table_type>)
      .def("insert", &insert<column_table_type>)
      .def("extend", &extend<column_table_type>)
      .def("update", &update<column_table_type>)
      .def("reorder", &reorder_table<column_table_type>)
      .def("sort", &sort_table<column_table_type>, (arg("column"), arg("reverse")=false))
      .def("types", &types<column_table_type>)
      ;

    boost::mpl::for_each<typename column_types::types>(
      column_table_set_data_wrapper<column_table_type>(
        column_table_class));
  }

}}}} // namespace dials::framework::boost_python::column_table_suite
