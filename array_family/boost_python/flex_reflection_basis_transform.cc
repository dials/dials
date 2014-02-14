/*
 * flex_reflection_basis_transform.cc
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
#include <cmath>
#include <scitbx/vec3.h>
#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <scitbx/array_family/ref_reductions.h>
#include <scitbx/array_family/boost_python/ref_pickle_double_buffered.h>
#include <scitbx/array_family/boost_python/flex_pickle_double_buffered.h>
#include <dials/algorithms/reflection_basis/transform.h>

namespace dials { namespace af { namespace boost_python {

  using namespace boost::python;
  using namespace scitbx::af::boost_python;
  using algorithms::reflection_basis::transform::Forward;
  using algorithms::reflection_basis::transform::TransformSpec;
  using model::Shoebox;
  using scitbx::vec3;

  /**
   * Perform a list of transforms and construct a flex array
   */
  template <typename FloatType>
  typename af::flex< Forward<FloatType> >::type* from_shoeboxes(
      const TransformSpec<FloatType> &spec,
      const af::const_ref< vec3<double> > &s1,
      const af::const_ref< double > &phi,
      const af::const_ref< Shoebox<FloatType> > &shoebox) {
    DIALS_ASSERT(s1.size() == phi.size());
    DIALS_ASSERT(s1.size() == shoebox.size());
    af::shared< Forward<FloatType> > result(s1.size());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = Forward<FloatType>(spec, s1[i], phi[i], shoebox[i]);
    }
    return new typename af::flex< Forward<FloatType> >::type(
        result, result.size());
  }

  /**
   * A class to convert the transform class to a string for pickling
   */
  template <typename FloatType>
  struct reflection_basis_transform_to_string : pickle_double_buffered::to_string
  {
    using pickle_double_buffered::to_string::operator<<;

    typedef Forward<FloatType> forward_type;

    /** Initialise with the version for checking */
    reflection_basis_transform_to_string() {
      unsigned int version = 1;
      *this << version;
    }

    /** Convert a single shoebox instance to string */
    reflection_basis_transform_to_string& operator<<(const forward_type &val) {
      DIALS_ERROR("Not Implemented");
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
  template <typename FloatType>
  struct reflection_basis_transform_from_string : pickle_double_buffered::from_string
  {
    using pickle_double_buffered::from_string::operator>>;

    typedef Forward<FloatType> forward_type;

    /** Initialise the class with the string. Get the version and check */
    reflection_basis_transform_from_string(const char* str_ptr)
    : pickle_double_buffered::from_string(str_ptr) {
      *this >> version;
      DIALS_ASSERT(version == 1);
    }

    /** Get a single shoebox instance from a string */
    reflection_basis_transform_from_string& operator>>(forward_type &val) {
      DIALS_ERROR("Not Implemented");
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

  template <typename FloatType>
  typename scitbx::af::boost_python::flex_wrapper<
    Forward<FloatType>,
    return_internal_reference<> >::class_f_t
  flex_reflection_basis_transform_wrapper(const char *name)
  {
    typedef Forward<FloatType> forward_type;

    return scitbx::af::boost_python::flex_wrapper <
      forward_type, return_internal_reference<> >::plain(name)
        .def("__init__", make_constructor(
          from_shoeboxes<FloatType>,
          default_call_policies(), (
            boost::python::arg("spec"),
            boost::python::arg("s1"),
            boost::python::arg("phi"),
            boost::python::arg("shoebox"))))
        .def_pickle(flex_pickle_double_buffered<forward_type,
          reflection_basis_transform_to_string<FloatType>,
          reflection_basis_transform_from_string<FloatType> >())
        ;
  }

  void export_flex_reflection_basis_transform() {
    typedef Shoebox<>::float_type float_type;
    flex_reflection_basis_transform_wrapper<float_type>(
        "reflection_basis_transform");
  }

}}} // namespace dials::af::boost_python
