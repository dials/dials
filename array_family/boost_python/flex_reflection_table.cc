/*
 * flex_reflection_table.cc
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
#include <dials/array_family/boost_python/flex_table_suite.h>
#include <dials/model/data/shoebox.h>
#include <dials/model/data/observation.h>
#include <scitbx/array_family/tiny_types.h>
#include <scitbx/array_family/ref_reductions.h>
#include <scitbx/vec3.h>
#include <scitbx/vec2.h>
#include <cctbx/miller.h>

namespace dials { namespace af { namespace boost_python {

  using namespace boost::python;
  using scitbx::vec2;
  using scitbx::vec3;
  using scitbx::af::int6;
  using flex_table_suite::flex_table_wrapper;
  using dials::model::Shoebox;
  using dials::model::Observation;

  /**
   * Construct a reflection table from a list of observations and shoeboxes
   * @param o The observation
   * @param s The shoeboxes
   * @returns A reflection table
   */
  template <typename T>
  T* make_from_observation_and_shoebox(
      const af::const_ref<Observation> &o,
      const af::const_ref< Shoebox<double> > &s) {
    DIALS_ASSERT(o.size() == s.size());

    // The reflection table
    T result(o.size());
    af::shared<std::size_t>    panel  = result["panel"];
    af::shared< vec3<double> > xyzval = result["xyzobs.px.value"];
    af::shared< vec3<double> > xyzvar = result["xyzobs.px.variance"];
    af::shared<double>         iraw   = result["intensity.raw.value"];
    af::shared<double>         irawv  = result["intensity.raw.variance"];
    af::shared<double>         icor   = result["intensity.cor.value"];
    af::shared<double>         icorv  = result["intensity.cor.variance"];
    af::shared< Shoebox<> >    sbox   = result["shoebox"];

    // Copy all the values
    for (std::size_t i = 0; i < result.nrows(); ++i) {

      // Check panel numbers
      DIALS_ASSERT(o[i].panel == s[i].panel);
      panel[i] = o[i].panel;

      // Copy observation info
      xyzval[i] = o[i].centroid.px.position;
      xyzvar[i] = o[i].centroid.px.std_err_sq;
      iraw[i]   = o[i].intensity.observed.value;
      irawv[i]  = o[i].intensity.observed.variance;
      icor[i]   = o[i].intensity.corrected.value;
      icorv[i]  = o[i].intensity.corrected.variance;

      // Copy shoebox info
      sbox[i].bbox = s[i].bbox;
      sbox[i].data = s[i].data;
      sbox[i].mask = s[i].mask;
      sbox[i].background = s[i].background;
    }

    // Return the new reflection table
    return new T(result);
  }


  /**
   * Struct to facilitate wrapping reflection table type
   */
  template <typename T>
  struct flex_reflection_table_wrapper : public flex_table_wrapper<T> {

    typedef flex_table_wrapper<T> base_type;
    typedef typename base_type::flex_types flex_types;
    typedef typename base_type::flex_table_type flex_table_type;
    typedef typename base_type::class_type class_type;

    /**
     * Wrap the reflection table class
     */
    static
    class_type wrap(const char *name) {

      // Wrap with flex table bindings
      class_type result = base_type::wrap(name);

      // Add functions
      result
        .def("__init__", make_constructor(
          &make_from_observation_and_shoebox<flex_table_type>))
        ;

      // return the wrapped class
      return result;
    }
  };

  void export_flex_reflection_table() {

    // Define all the types we want to support in the table
    typedef flex_type_generator<
      bool,
      int,
      std::size_t,
      double,
      std::string,
      vec2<double>,
      vec3<double>,
      int6,
      cctbx::miller::index<>,
      Shoebox<>
    >::type flex_types;

    // Export the reflection table
    flex_reflection_table_wrapper<flex_types>::wrap("reflection_table");
  }

}}} // namespace dials::af::boost_python
