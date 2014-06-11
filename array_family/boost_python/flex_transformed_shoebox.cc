/*
 * flex_transformed_shoebox.cc
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
#include <dials/algorithms/reflection_basis/ideal_profile.h>
#include <dials/model/data/transformed_shoebox.h>

namespace dials { namespace af { namespace boost_python {

  using namespace boost::python;
  using namespace scitbx::af::boost_python;
  using algorithms::reflection_basis::transform::Forward;
  using algorithms::reflection_basis::transform::TransformSpec;
  using algorithms::reflection_basis::ideal_profile;
  using model::Shoebox;
  using model::TransformedShoebox;
  using scitbx::vec3;

  /**
   * Perform a list of transforms and construct a flex array
   */
  static
  af::flex<TransformedShoebox>::type* from_shoeboxes(
      const TransformSpec<> &spec,
      const af::const_ref< vec3<double> > &s1,
      const af::const_ref< double > &phi,
      const af::const_ref< Shoebox<> > &shoebox) {
    DIALS_ASSERT(s1.size() == phi.size());
    DIALS_ASSERT(s1.size() == shoebox.size());
    af::shared<TransformedShoebox> result(s1.size());
    for (std::size_t i = 0; i < result.size(); ++i) {
      Forward<> transform(spec, s1[i], phi[i], shoebox[i]);
      result[i].data = transform.profile();
      result[i].background = transform.background();
    }
    return new af::flex<TransformedShoebox>::type(result, result.size());
  }

  /**
   * Compute the correlation between profiles
   */
  template <typename FloatType>
  FloatType
  compute_correlation(const af::const_ref<FloatType, af::c_grid<3> > &p,
                      const af::const_ref<FloatType, af::c_grid<3> > &c,
                      const af::const_ref<FloatType, af::c_grid<3> > &b) {
    double xb = 0.0, yb = 0.0;
    std::size_t count = 0;
    for (std::size_t i = 0; i < p.size(); ++i) {
      xb += p[i];
      yb += c[i] - b[i];
      count++;
    }
    DIALS_ASSERT(count > 0);
    xb /= count;
    yb /= count;
    double sdxdy = 0.0, sdx2 = 0.0, sdy2 = 0.0;
    for (std::size_t i = 0; i < p.size(); ++i) {
      double dx = p[i] - xb;
      double dy = c[i] - b[i] - yb;
      sdxdy += dx*dy;
      sdx2 += dx*dx;
      sdy2 += dy*dy;
    }
    if (sdx2 == 0.0 || sdy2 == 0.0) {
      return 0.0;
    }
    return sdxdy / (std::sqrt(sdx2) * std::sqrt(sdy2));
  }

  /**
   * Compute the correlation of the profiles with the idea profile
   */
  static
  af::shared<double> correlation_with_ideal(
      const af::const_ref<TransformedShoebox> &self,
      double n_sigma) {
    af::shared<double> result(self.size());

    // Get the size of the profile
    DIALS_ASSERT(self.size() > 0);
    std::size_t size = self[0].data.accessor()[0] / 2;

    // Create the ideal profile
    af::versa< double, af::c_grid<3> > ideal = ideal_profile<double>(size, n_sigma);

    // Compute the correlations
    for (std::size_t i = 0; i < self.size(); ++i) {
      DIALS_ASSERT(ideal.accessor().all_eq(self[i].data.accessor()));
      DIALS_ASSERT(ideal.accessor().all_eq(self[i].background.accessor()));
      result[i] = compute_correlation(
          ideal.const_ref(),
          self[i].data.const_ref(),
          self[i].background.const_ref());
    }

    // Return the array
    return result;
  }

  /**
   * A class to convert the transform class to a string for pickling
   */
  struct transformed_shoebox_to_string : pickle_double_buffered::to_string
  {
    using pickle_double_buffered::to_string::operator<<;

    /** Initialise with the version for checking */
    transformed_shoebox_to_string() {
      unsigned int version = 1;
      *this << version;
    }

    /** Convert a single shoebox instance to string */
    transformed_shoebox_to_string& operator<<(const TransformedShoebox &val) {
      profile_to_string(val.data);
      profile_to_string(val.background);
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
  struct transformed_shoebox_from_string : pickle_double_buffered::from_string
  {
    using pickle_double_buffered::from_string::operator>>;

    /** Initialise the class with the string. Get the version and check */
    transformed_shoebox_from_string(const char* str_ptr)
    : pickle_double_buffered::from_string(str_ptr) {
      *this >> version;
      DIALS_ASSERT(version == 1);
    }

    /** Get a single shoebox instance from a string */
    transformed_shoebox_from_string& operator>>(TransformedShoebox &val) {
      typedef af::versa< double, af::c_grid<3> > profile_type;
      val.data = profile_from_string<profile_type>();
      val.background = profile_from_string<profile_type>();
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
    TransformedShoebox,
    return_internal_reference<> >::class_f_t
  flex_transformed_shoebox_wrapper(const char *name)
  {
    return scitbx::af::boost_python::flex_wrapper <
      TransformedShoebox, return_internal_reference<> >::plain(name)
        .def("__init__", make_constructor(
          from_shoeboxes,
          default_call_policies(), (
            boost::python::arg("spec"),
            boost::python::arg("s1"),
            boost::python::arg("phi"),
            boost::python::arg("shoebox"))))
        .def("correlation_with_ideal",
            &correlation_with_ideal)
        .def_pickle(flex_pickle_double_buffered<
          TransformedShoebox,
          transformed_shoebox_to_string,
          transformed_shoebox_from_string>())
        ;
  }

  void export_flex_transformed_shoebox() {
    flex_transformed_shoebox_wrapper("transformed_shoebox");
  }

}}} // namespace dials::af::boost_python
