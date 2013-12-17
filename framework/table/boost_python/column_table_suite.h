/*
 * column_table_suite.h
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





  /**
   * A visitor to convert the column to a boost python object
   */
  struct column_to_object_visitor : public boost::static_visitor<object> {
    template <typename T>
    object operator () (T &col) {
      return object(col);
    }
  };

  /**
   * A visitor to set the row item from a python object
   */
  struct setitem_row_visitor : public boost::static_visitor<void> {
    std::size_t index;
    object item;
    
    setitem_row_visitor(std::size_t index_, object item_)
      : index(index_),
        item(item_) {}

    template <typename T>
    void operator () (T &col) {
      col[index] = extract<typename T::value_type>(item);
    }
  };

  /**
   * A visitor to append column data from 1 table to another
   */
  template <typename T>
  struct append_column_visitor : public boost::static_visitor<void> {

    T &self;
    typename T::key_type key;
    typename T::size_type na, nb;

    append_column_visitor(
          T &self_,
          typename T::key_type key_,
          typename T::size_type na_, 
          typename T::size_type nb_)
      : self(self_),
        key(key_),
        na(na_),
        nb(nb_) {}

    template <typename U>
    void operator () (const U &col) {
      U c = self[key];
      for (typename T::size_type i = 0; i < nb; ++i) {
        c[na + i] = col[i];
      }
    }
  };
  
  /**
   * A visitor to add new columns (and over-write old columns) in the table.
   */
  template <typename T>
  struct update_column_visitor : public boost::static_visitor<void> {

    T &self;
    typename T::key_type key;

    update_column_visitor(T &self_, typename T::key_type key_)
      : self(self_),
        key(key_) {}

    template <typename U>
    void operator () (const U &col) {
      self.erase(key);
      U c = self[key];
      for (std::size_t i = 0; i < col.size(); ++i) {
        c[i] = col[i];
      }
    }
  };
  
  /**
   * Read a slice from the input column and write it to a column in the
   * output table.
   */
  template <typename T>
  struct copy_to_slice_visitor : public boost::static_visitor<void> {

    T &self;
    typename T::key_type key;
    scitbx::boost_python::adapted_slice slice;

    copy_to_slice_visitor(
          T &self_,
          typename T::key_type key_,
          scitbx::boost_python::adapted_slice slice_)
      : self(self_),
        key(key_),
        slice(slice_) {}

    template <typename U>
    void operator () (const U &col) {
      U c = self[key];
      for (std::size_t i = 0, j = slice.start; 
          i < self.nrows(); ++i, j += slice.step) {
        c[i] = col[j];
      }
    }
  };
  
  /**
   * Copy from a slice of column data into the column table.
   */
  template <typename T>
  struct copy_from_slice_visitor : public boost::static_visitor<void> {

    T &self;
    typename T::key_type key;
    scitbx::boost_python::adapted_slice slice;
    typename T::size_type num;
    
    copy_from_slice_visitor(
          T &self_,
          typename T::key_type key_,
          scitbx::boost_python::adapted_slice slice_,
          std::size_t num_)
      : self(self_),
        key(key_),
        slice(slice_),
        num(num_) {}

    template <typename U>
    void operator () (const U &col) {
      U c = self[key];
      for (std::size_t i = 0, j = slice.start; 
          i < num; ++i, j += slice.step) {
        c[j] = col[i];
      }
    }
  };
  
 
  /**
   * A visitor to reorder the elements of a column
   */
  struct reorder_visitor : public boost::static_visitor<void> {

    af::const_ref<std::size_t> index;

    reorder_visitor(const af::const_ref<std::size_t> &index_)
      : index(index_) {}

    template <typename T>
    void operator () (T &col) {
      std::vector<typename T::value_type> temp(col.begin(), col.end());
      for (std::size_t i = 0; i < index.size(); ++i) {
        col[i] = temp[index[i]];
      }
    }
  };
  
  /**
   * Functor to compare elements by index
   */
  template <typename T>
  struct compare_index {
    const T &v_;

    compare_index(const T &v)
      : v_(v) {}
      
    template <typename U>
    bool operator()(U a, U b) {
      return v_[a] < v_[b];
    }
  };

  /**
   * A visitor to sort the table by columns
   */
  struct sort_visitor : public boost::static_visitor<void> {
    af::ref<std::size_t> index;
    
    sort_visitor(af::ref<std::size_t> index_)
        : index(index_) {
      for (std::size_t i = 0; i < index.size(); ++i) {
        index[i] = i;
      }
    }

    template <typename T>
    void operator () (const T &col) {
      std::sort(index.begin(), index.end(), compare_index<T>(col));
    }
  };
  
  /**
   * Initialise the column table from a list of (key, column) pairs
   * @param columns The list of columns
   * @returns The column table
   */
  template <typename T>
  T* make_column_table(object columns) {
    T *self = new T();
    object obj(self);
    for (std::size_t i = 0; i < len(columns); ++i) {
      obj[columns[i][0]] = columns[i][1];
    }
    return self;
  }

  /**
   * Get a column of data
   * @param self The column table
   * @param key The column name
   * @returns The boost python object with the column data
   */
  template <typename T>
  object getitem_column(T &self, const typename T::key_type &key) {
    typename T::mapped_type col = self[key];
    DIALS_ASSERT(!col.empty());
    column_to_object_visitor visitor;
    return col.apply_visitor(visitor);
  }

  /**
   * Set a column of data
   * @param self The table to populate
   * @param key The column name
   * @param data The column data
   */
  template <typename T, typename U>
  void setitem_column(
      T &self,
      const typename T::key_type &key,
      const af::const_ref<U> &data) {
    self.erase(key);
    column_data<U> col = self[key];
    col.resize(data.size());
    std::copy(data.begin(), data.end(), col.begin());
  }

  /**
   * An MPL generator to create setitem functions for each type
   */
  template <typename T>
  struct setitem_column_generator {
    class_<T> table_class;
    setitem_column_generator(class_<T> table_class_)
      : table_class(table_class_) {}
    template <typename U>
    void operator()(const U &x) {
      table_class.def("__setitem__",
        &setitem_column<T, typename U::value_type>,
          return_internal_reference<>());
    }
  };

  /**
   * Get a row of data from the table.
   * @param self The table
   * @param n The position of the row
   */
  template <typename T>
  dict getitem_row(const T &self, typename T::size_type n) {
    typedef typename T::const_iterator iterator; 
    dict result;
    column_to_object_visitor visitor;
    for (iterator it = self.begin(); it != self.end(); ++it) {
      result[it->first] = it->second.apply_visitor(visitor)[n];
    }
    return result;
  }

  /**
   * Set a row of data in the table
   * @param self The table to modify
   * @param n The position of the row
   * @param row The row data to set
   */
  template <typename T>
  void setitem_row(T &self, typename T::size_type n, dict row) {
    for (typename T::iterator it = self.begin(); it != self.end(); ++it) {
      if (row.has_key(it->first)) {
        setitem_row_visitor visitor(n, row[it->first]);
        it->second.apply_visitor(visitor);
      }
    }
  }

  /**
   * Get a slice of the table and return a new table
   * @param self The current table
   * @param slice The slice
   * @returns A new table with the chosen elements
   */
  template <typename T>
  T getitem_slice(const T &self, slice s) {
    typedef typename T::const_iterator iterator; 
    scitbx::boost_python::adapted_slice as(s, self.nrows());
    T result(as.size);
    for (iterator it = self.begin(); it != self.end(); ++it) {
      copy_to_slice_visitor<T> visitor(result, it->first, as);
      it->second.apply_visitor(visitor);
    }
    return result;
  }

  /**
   * Set a slice of data from one table into the current table.
   * @param self The current table
   * @param slice The slice to set
   * @param other The other table whose elements to set
   */
  template <typename T>
  void setitem_slice(T &self, slice s, const T &other) {
    typedef typename T::const_iterator iterator; 
    scitbx::boost_python::adapted_slice as(s, self.nrows());
    for (iterator it = other.begin(); it != other.end(); ++it) {
      copy_from_slice_visitor<T> visitor(self, it->first, as, other.nrows());
      it->second.apply_visitor(visitor);
    }
  }

  /**
   * Check if the table has the given key
   * @param self The table
   * @param key The key to check for
   */
  template <typename T>
  bool has_key(const T &self, const typename T::key_type &key) {
    return self.count(key) == 1;
  }

  /**
   * Append a row onto the end of the columns
   * @param self The column table
   * @param row The row to add
   */
  template <typename T>
  void append(T &self, dict row) {
    self.resize(self.nrows() + 1);
    setitem_row(self, self.nrows() - 1, row);
  }

  /**
   * Insert a row into the table at the given position
   * @param self The column table
   * @param n The position to insert the row at
   * @param row The row to insert
   */
  template <typename T>
  void insert(T &self, typename T::size_type n, dict row) {
    self.insert(n);
    setitem_row(self, n, row);
  }

  /**
   * Extend the table with column data from another table. This will add
   * all the columns from the other table onto the end of the current table.
   * @param self The current table
   * @param other The other table
   */
  template <typename T>
  void extend(T &self, const T &other) {
    typedef typename T::const_iterator iterator; 
    typename T::size_type ns = self.nrows();
    typename T::size_type no = other.nrows();
    self.resize(ns + no);
    for (iterator it = other.begin(); it != other.end(); ++it) {
      append_column_visitor<T> visitor(self, it->first, ns, no);
      it->second.apply_visitor(visitor);
    }
  }

  /**
   * Update the table with column data from another table. New columns are added
   * to the table and exisiting columns are over-written by columns from the
   * other table.
   * @param self The current table
   * @param other The other table
   */
  template <typename T>
  void update(T &self, const T &other) {
    typedef typename T::const_iterator iterator; 
    for (iterator it = other.begin(); it != other.end(); ++it) {
      update_column_visitor<T> visitor(self, it->first);
      it->second.apply_visitor(visitor);
    }
  }

  /**
   * An MPL function to populate the list of types
   */
  struct type_appender {
    list type_list;
    
    type_appender(list type_list_)
      : type_list(type_list_) {}
    
    template <typename U>
    void operator()(U x) {
      typename U::value_type a;
      type_list.append(object(handle<>(PyObject_Type(object(a).ptr()))));
    }
  };

  /**
   * Get a list of valid column types.
   * @param self The table
   * @returns A list of python types.
   */
  template <typename T>
  list types(const T &self) {
    list result;
    boost::mpl::for_each<typename T::mapped_type::types>(type_appender(result));
    return result;
  }

  /**
   * Reorder all the columns according to the input indices
   * @param self The table object
   * @param indices The array of indices
   */
  template <typename T>
  void reorder(T &self, const af::const_ref<std::size_t> &index) {
    typedef typename T::iterator iterator; 
    reorder_visitor visitor(index);
    for (iterator it = self.begin(); it != self.end(); ++it) {
      it->second.apply_visitor(visitor);
    }
  }

  /**
   * Sort the table with respect to a given column
   * @param self The table object
   * @param key The column key
   * @param reverse True/False reverse the sort
   */
  template <typename T>
  void sort(T &self, typename T::key_type key, bool reverse) {
    af::shared<std::size_t> index(self.nrows());
    sort_visitor visitor(index.ref());
    typename T::mapped_type col = self[key];
    col.apply_visitor(visitor);
    if (reverse) {
      std::reverse(index.begin(), index.end());
    }
    reorder(self, index.const_ref());
  }
  
  
  /**
   * A proxy iterator to iterate over the table keys
   */
  template <typename T>
  class key_iterator {
  public:

    typedef std::forward_iterator_tag iterator_category;
    typedef T table_type;
    typedef typename T::const_iterator base_iterator;
    typedef ptrdiff_t difference_type;
    typedef typename T::key_type value_type;
    typedef const value_type *pointer;
    typedef const value_type &reference;

    key_iterator(base_iterator it)
      : it_(it) {}

    reference operator*() {
      return it_->first;
    }

    key_iterator& operator++() {
      ++it_;
      return *this;
    }

    key_iterator operator++(int) {
      key_iterator result(*this);
      ++(*this);
      return result;
    }

    bool operator==(const key_iterator& rhs) const {
      return it_ == rhs.it_;
    }

    bool operator!=(const key_iterator& rhs) const {
      return !(*this == rhs);
    }

  private:
    base_iterator it_;
  };


  /**
   * A proxy iterator to iterate over the table column
   */
  template <typename T>
  class column_iterator {
  public:

    typedef T table_type;
    typedef typename T::const_iterator base_iterator;
    typedef ptrdiff_t difference_type;
    typedef std::forward_iterator_tag iterator_category;
    typedef object value_type;
    typedef const value_type *pointer;
    typedef const value_type reference;

    column_iterator(base_iterator it)
      : it_(it) {}

    reference operator*() {
      column_to_object_visitor visitor;
      return boost::python::make_tuple(
        it_->first,
        it_->second.apply_visitor(visitor));
    }

    column_iterator& operator++() {
      ++it_;
      return *this;
    }

    column_iterator operator++(int) {
      column_iterator result(*this);
      ++(*this);
      return result;
    }

    bool operator==(const column_iterator& rhs) const {
      return it_ == rhs.it_;
    }

    bool operator!=(const column_iterator& rhs) const {
      return !(*this == rhs);
    }

  private:
    base_iterator it_;
  };


  /**
   * Proxy iterator to iterate through the table rows
   */
  template <typename T>
  class row_iterator {
  public:

    typedef T table_type;
    typedef ptrdiff_t difference_type;
    typedef std::forward_iterator_tag iterator_category;
    typedef dict value_type;
    typedef const value_type *pointer;
    typedef const value_type reference;

    row_iterator(const T &self, std::size_t index)
      : self_(self),
        index_(index) {}

    reference operator*() {
      typedef typename T::const_iterator iterator; 
      dict result;
      column_to_object_visitor visitor;
      for (iterator it = self_.begin(); it != self_.end(); ++it) {
        result[it->first] = it->second.apply_visitor(visitor)[index_];
      }
      return result;
    }

    row_iterator& operator++() {
      ++index_;
      return *this;
    }

    row_iterator operator++(int) {
      row_iterator result(*this);
      ++(*this);
      return result;
    }

    bool operator==(const row_iterator& rhs) const {
      return index_ == rhs.index_;
    }

    bool operator!=(const row_iterator& rhs) const {
      return !(*this == rhs);
    }

  private:
    const T &self_;
    std::size_t index_;
  };


  /**
   * Struct to help in creation of table proxy iterators
   */
  template <typename Iterator>
  struct make_iterator {
    static
    Iterator begin(const typename Iterator::table_type &self) {
      return Iterator(self.begin());
    }
    
    static
    Iterator end(const typename Iterator::table_type &self) {
      return Iterator(self.end());
    }
  };

  /**
   * Specialization for row iterator
   */
  template <typename T>
  struct make_iterator< row_iterator<T> > {
    static
    row_iterator<T> begin(const T &self) {
      return row_iterator<T>(self, 0);
    }
    
    static
    row_iterator<T> end(const T &self) {
      return row_iterator<T>(self, self.nrows());
    }
  };
  
  /**
   * Create an iterator range for the table iterators
   */
  template <typename Iterator>
  object make_iterator_range() {
    return range(&make_iterator<Iterator>::begin,
                 &make_iterator<Iterator>::end);
  }
  
  
  /**
   * Export the wrapped class to python
   */
  template <typename column_types>
  void column_table_wrapper(const char *name) {

    boost::mpl::for_each<typename column_types::types>(
      column_data_wrapper(name));

    typedef column_table<column_types> column_table_type;

    class_<column_table_type> column_table_class(name);
    column_table_class
      .def(init<std::size_t>())
      .def("__init__", make_constructor(
        &make_column_table<column_table_type>))
      .def("types", &types<column_table_type>)
      .def("has_key", &has_key<column_table_type>)
      .def("clear", &column_table_type::clear)
      .def("empty", &column_table_type::empty)
      .def("append", &append<column_table_type>)
      .def("insert", &insert<column_table_type>)
      .def("extend", &extend<column_table_type>)
      .def("update", &update<column_table_type>)
      .def("nrows", &column_table_type::nrows)
      .def("ncols", &column_table_type::ncols)
      .def("__len__", &column_table_type::size)
      .def("__contains__", &has_key<column_table_type>)
      .def("__getitem__", &getitem_column<column_table_type>)
      .def("__getitem__", &getitem_row<column_table_type>)
      .def("__setitem__", &setitem_row<column_table_type>)
      .def("__getitem__", &getitem_slice<column_table_type>)
      .def("__setitem__", &setitem_slice<column_table_type>)
      .def("cols", make_iterator_range< column_iterator<column_table_type> >())
      .def("rows", make_iterator_range< row_iterator<column_table_type> >())
      .def("keys", make_iterator_range< key_iterator<column_table_type> >())
      .def("reorder", &reorder<column_table_type>)
      .def("sort", &sort<column_table_type>, (arg("column"), arg("reverse")=false))
      ;

    // For each column type, create a __setitem__ method to set column data
    boost::mpl::for_each<typename column_types::types>(
      setitem_column_generator<column_table_type>(column_table_class));
  }

}}}} // namespace dials::framework::boost_python::column_table_suite
