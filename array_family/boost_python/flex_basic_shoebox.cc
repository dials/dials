/*
 * flex_shoebox.cc
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
#include <boost/unordered_map.hpp>
#include <boost_adaptbx/std_pair_conversion.h>
#include <cmath>
#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <scitbx/array_family/ref_reductions.h>
#include <scitbx/array_family/boost_python/ref_pickle_double_buffered.h>
#include <scitbx/array_family/boost_python/flex_pickle_double_buffered.h>
#include <dials/model/data/basic_shoebox.h>
#include <dials/model/data/shoebox.h>
#include <dials/config.h>

namespace dials { namespace af { namespace boost_python {

  using namespace boost::python;
  using namespace scitbx::af::boost_python;

  using af::int2;
  using af::int6;
  using af::small;
  using dials::model::BasicShoebox;
  using dials::model::Shoebox;

  /**
   * Construct from an array of panels and bounding boxes.
   */
  static
  af::flex<BasicShoebox>::type* from_panel_and_bbox(
      const af::const_ref<std::size_t> panel,
      const af::const_ref<int6> bbox) {
    DIALS_ASSERT(panel.size() == bbox.size());
    af::shared<BasicShoebox> result(panel.size());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = BasicShoebox(panel[i], bbox[i]);
    }
    return new af::flex<BasicShoebox>::type(
      result, af::flex_grid<>(result.size()));
  }

  /**
   * Check if the arrays are consistent
   */
  static
  shared<bool> is_consistent(const const_ref<BasicShoebox> &a) {
    shared<bool> result(a.size(), af::init_functor_null<bool>());
    for (std::size_t i = 0; i < a.size(); ++i) {
      result[i] = a[i].is_consistent();
    }
    return result;
  }

  /**
   * Get the bounding boxes
   */
  static
  shared<int6> bounding_boxes(const const_ref<BasicShoebox> &a) {
    shared<int6> result(a.size(), af::init_functor_null<int6>());
    for (std::size_t i = 0; i < a.size(); ++i) {
      result[i] = a[i].bbox;
    }
    return result;
  }

  /**
   * Get the panel numbers
   */
  static
  shared<std::size_t> panels(const const_ref<BasicShoebox> &a) {
    shared<std::size_t> result(a.size(), af::init_functor_null<std::size_t>());
    for (std::size_t i = 0; i < a.size(); ++i) {
      result[i] = a[i].panel;
    }
    return result;
  }

  /**
   * A class to convert the shoebox class to a string for pickling
   */
  struct basic_shoebox_to_string : pickle_double_buffered::to_string
  {
    using pickle_double_buffered::to_string::operator<<;

    /** Initialise with the version for checking */
    basic_shoebox_to_string() {
      unsigned int version = 1;
      *this << version;
    }

    /** Convert a single shoebox instance to string */
    basic_shoebox_to_string& operator<<(const BasicShoebox &val) {
      *this << val.panel
            << val.bbox[0]
            << val.bbox[1]
            << val.bbox[2]
            << val.bbox[3]
            << val.bbox[4]
            << val.bbox[5];

      profile_to_string(val.data);

      return *this;
    }

    /** Convert a profile to string */
    template <typename ProfileType>
    void profile_to_string(const ProfileType &p) {
      *this << p.accessor().size();
      for (std::size_t i = 0; i < p.accessor().size(); ++i) {
        *this << p.accessor()[i];
      }
      for (std::size_t i = 0; i < p.size(); ++i) {
        *this << p[i];
      }
    }
  };

  /**
   * A class to convert a string to a shoebox for unpickling
   */
  struct basic_shoebox_from_string : pickle_double_buffered::from_string
  {
    using pickle_double_buffered::from_string::operator>>;

    /** Initialise the class with the string. Get the version and check */
    basic_shoebox_from_string(const char* str_ptr)
    : pickle_double_buffered::from_string(str_ptr) {
      *this >> version;
      DIALS_ASSERT(version == 1);
    }

    /** Get a single shoebox instance from a string */
    basic_shoebox_from_string& operator>>(BasicShoebox &val) {
      *this >> val.panel
            >> val.bbox[0]
            >> val.bbox[1]
            >> val.bbox[2]
            >> val.bbox[3]
            >> val.bbox[4]
            >> val.bbox[5];

      val.data = profile_from_string< versa<int, c_grid<3> > >();

      return *this;
    }

    /** Get a profile from a string */
    template <typename ProfileType>
    ProfileType profile_from_string() {
      typename ProfileType::accessor_type accessor;
      typename ProfileType::size_type n_dim;
      *this >> n_dim;
      DIALS_ASSERT(n_dim == accessor.size());
      for (std::size_t i = 0; i < n_dim; ++i) {
        *this >> accessor[i];
      }
      ProfileType p = ProfileType(accessor);
      for (std::size_t i = 0; i < p.size(); ++i) {
        *this >> p[i];
      }
      return p;
    }

    unsigned int version;
  };

  scitbx::af::boost_python::flex_wrapper<
    BasicShoebox,
    return_internal_reference<> >::class_f_t
  flex_basic_shoebox_wrapper(const char *name)
  {
    return scitbx::af::boost_python::flex_wrapper <
      BasicShoebox, return_internal_reference<> >::plain(name)
        .def("__init__", make_constructor(
          from_panel_and_bbox,
          default_call_policies(), (
            boost::python::arg("panel"),
            boost::python::arg("bbox"))))
        .def("is_consistent", &is_consistent)
        .def("panels", &panels)
        .def("bounding_boxes", &bounding_boxes)
        .def_pickle(flex_pickle_double_buffered<BasicShoebox,
          basic_shoebox_to_string,
          basic_shoebox_from_string>());
  }

  void export_flex_basic_shoebox() {
    flex_basic_shoebox_wrapper("basic_shoebox");

    boost_adaptbx::std_pair_conversions::to_and_from_tuple<
      af::shared<std::size_t>,
      af::shared<BasicShoebox> >();
  }

}}} // namespace dials::af::boost_python
