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
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/mpl/for_each.hpp>
#include <string>
#include <scitbx/array_family/flex_types.h>
#include <dials/framework/table/column_table.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/error.h>

namespace dials { namespace framework { namespace boost_python {

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
      table_class_.def("__setitem__", &column_table_set_data<T, typename U::value_type>,
        return_internal_reference<>());
    }
  };

  struct column_data_to_object : public boost::static_visitor<boost::python::object> {

    template <typename T>
    boost::python::object operator () (T &col)
    {
      return boost::python::object(col);
    }
  };

  template <typename T>
  boost::python::object column_table_get_data(
      column_table<T> &table,
      const typename column_table<T>::key_type &key) {

    typedef typename column_table<T>::mapped_type mapped_type;
    mapped_type col = table[key];
    DIALS_ASSERT(!col.empty());
    column_data_to_object to_object;
    return col.apply_visitor(to_object);
  }


  template <typename column_types>
  void column_table_wrapper(const char *name) {

    boost::mpl::for_each<typename column_types::types>(
      column_data_wrapper(name));


    typedef column_table<column_types> column_table_type;
    class_<column_table_type> column_table_class(name);
    column_table_class
      .def("erase", &column_table_type::erase)
      .def("empty", &column_table_type::empty)
      .def("clear", &column_table_type::clear)
      .def("nrows", &column_table_type::nrows)
      .def("ncols", &column_table_type::ncols)
      .def("__len__", &column_table_type::size)
      .def("__getitem__", &column_table_get_data<typename column_table_type::mapped_type>);

    boost::mpl::for_each<typename column_types::types>(
      column_table_set_data_wrapper<column_table_type>(column_table_class));
  }

  void export_column_table() {

    typedef column_type_generator<
      int,
      double,
      std::string
    >::type column_types;

    column_table_wrapper<column_types>("column_table");

  }

}}} // namespace dials::framework::boost_python
