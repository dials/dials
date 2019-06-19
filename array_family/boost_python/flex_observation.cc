/*
 * flex_observation.cc
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
#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <scitbx/array_family/ref_reductions.h>
#include <scitbx/array_family/boost_python/ref_pickle_double_buffered.h>
#include <scitbx/array_family/boost_python/flex_pickle_double_buffered.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/model/data/observation.h>
#include <dials/error.h>

namespace dials { namespace af { namespace boost_python {

  using namespace boost::python;
  using namespace scitbx::af::boost_python;

  using dials::model::Centroid;
  using dials::model::Intensity;
  using dials::model::Observation;
  using dxtbx::model::BeamBase;
  using dxtbx::model::Detector;
  using dxtbx::model::Scan;
  using scitbx::vec2;
  using scitbx::vec3;

  /** Initialise an observation list from panels */
  static af::flex<Observation>::type *init_from_panel(
    const af::const_ref<std::size_t> &panel) {
    af::shared<Observation> result(panel.size(), Observation());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i].panel = panel[i];
    }
    return new af::flex<Observation>::type(result, af::flex_grid<>(result.size()));
  }

  /** Initialise an observation list from centroids */
  static af::flex<Observation>::type *init_from_centroid(
    const af::const_ref<Centroid> &centroid) {
    af::shared<Observation> result(centroid.size(), Observation());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i].centroid = centroid[i];
    }
    return new af::flex<Observation>::type(result, af::flex_grid<>(result.size()));
  }

  /** Initialise an observation list from intensities */
  static af::flex<Observation>::type *init_from_intensity(
    const af::const_ref<Intensity> &intensity) {
    af::shared<Observation> result(intensity.size(), Observation());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i].intensity = intensity[i];
    }
    return new af::flex<Observation>::type(result, af::flex_grid<>(result.size()));
  }

  /** Initialise an observation list from centroids and intensities */
  static af::flex<Observation>::type *init_from_centroid_and_intensity(
    const af::const_ref<Centroid> &centroid,
    const af::const_ref<Intensity> &intensity) {
    DIALS_ASSERT(centroid.size() == intensity.size());
    af::shared<Observation> result(intensity.size(), Observation());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i].centroid = centroid[i];
      result[i].intensity = intensity[i];
    }
    return new af::flex<Observation>::type(result, af::flex_grid<>(result.size()));
  }

  /** Initialise an observation list from panels and centroids */
  static af::flex<Observation>::type *init_from_panel_and_centroid(
    const af::const_ref<std::size_t> &panel,
    const af::const_ref<Centroid> &centroid) {
    DIALS_ASSERT(panel.size() == centroid.size());
    af::shared<Observation> result(centroid.size(), Observation());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i].panel = panel[i];
      result[i].centroid = centroid[i];
    }
    return new af::flex<Observation>::type(result, af::flex_grid<>(result.size()));
  }

  /** Initialise an observation list from panels and intensities */
  static af::flex<Observation>::type *init_from_panel_and_intensity(
    const af::const_ref<std::size_t> &panel,
    const af::const_ref<Intensity> &intensity) {
    DIALS_ASSERT(panel.size() == intensity.size());
    af::shared<Observation> result(intensity.size(), Observation());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i].panel = panel[i];
      result[i].intensity = intensity[i];
    }
    return new af::flex<Observation>::type(result, af::flex_grid<>(result.size()));
  }

  /** Initialise an observation list from panels, centroids and intensities */
  static af::flex<Observation>::type *init_from_panel_centroid_and_intensity(
    const af::const_ref<std::size_t> &panel,
    const af::const_ref<Centroid> &centroid,
    const af::const_ref<Intensity> &intensity) {
    DIALS_ASSERT(centroid.size() == intensity.size());
    DIALS_ASSERT(panel.size() == intensity.size());
    af::shared<Observation> result(intensity.size(), Observation());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i].panel = panel[i];
      result[i].centroid = centroid[i];
      result[i].intensity = intensity[i];
    }
    return new af::flex<Observation>::type(result, af::flex_grid<>(result.size()));
  }

  /** Initialise an observation list from panel and centroids */
  static af::flex<Observation>::type *init_from_single_panel_and_centroid(
    std::size_t panel,
    const af::const_ref<Centroid> &centroid) {
    af::shared<Observation> result(centroid.size(), Observation());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i].panel = panel;
      result[i].centroid = centroid[i];
    }
    return new af::flex<Observation>::type(result, af::flex_grid<>(result.size()));
  }

  /** Initialise an observation list from panel and intensities */
  static af::flex<Observation>::type *init_from_single_panel_and_intensity(
    std::size_t panel,
    const af::const_ref<Intensity> &intensity) {
    af::shared<Observation> result(intensity.size(), Observation());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i].panel = panel;
      result[i].intensity = intensity[i];
    }
    return new af::flex<Observation>::type(result, af::flex_grid<>(result.size()));
  }

  /** Initialise an observation list from panel, centroids and intensities */
  static af::flex<Observation>::type *init_from_single_panel_centroid_and_intensity(
    std::size_t panel,
    const af::const_ref<Centroid> &centroid,
    const af::const_ref<Intensity> &intensity) {
    DIALS_ASSERT(centroid.size() == intensity.size());
    af::shared<Observation> result(intensity.size(), Observation());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i].panel = panel;
      result[i].centroid = centroid[i];
      result[i].intensity = intensity[i];
    }
    return new af::flex<Observation>::type(result, af::flex_grid<>(result.size()));
  }

  /** @returns An array of panel numbers */
  static af::shared<std::size_t> observation_get_panels(
    const af::const_ref<Observation> &obj) {
    af::shared<std::size_t> result(obj.size(), af::init_functor_null<std::size_t>());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = obj[i].panel;
    }
    return result;
  }

  /** @returns An array of centroids */
  static af::shared<Centroid> observation_get_centroids(
    const af::const_ref<Observation> &obj) {
    af::shared<Centroid> result(obj.size(), Centroid());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = obj[i].centroid;
    }
    return result;
  }

  /** @returns An array of intensities */
  static af::shared<Intensity> observation_get_intensities(
    const af::const_ref<Observation> &obj) {
    af::shared<Intensity> result(obj.size(), Intensity());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = obj[i].intensity;
    }
    return result;
  }

  /** @returns An array of panel numbers */
  static void observation_set_panels(af::ref<Observation> obj,
                                     const af::const_ref<std::size_t> &panel) {
    DIALS_ASSERT(obj.size() == panel.size());
    for (std::size_t i = 0; i < obj.size(); ++i) {
      obj[i].panel = panel[i];
    }
  }

  /** @returns An array of centroids */
  static void observation_set_centroids(af::ref<Observation> obj,
                                        const af::const_ref<Centroid> &centroid) {
    DIALS_ASSERT(obj.size() == centroid.size());
    for (std::size_t i = 0; i < obj.size(); ++i) {
      obj[i].centroid = centroid[i];
    }
  }

  /** @returns An array of intensities */
  static void observation_set_intensities(af::ref<Observation> obj,
                                          const af::const_ref<Intensity> &intensity) {
    DIALS_ASSERT(obj.size() == intensity.size());
    for (std::size_t i = 0; i < obj.size(); ++i) {
      obj[i].intensity = intensity[i];
    }
  }

  /** Update the millimeter centroid positions of all observations */
  void observation_update_centroid_mm(af::ref<Observation> obj,
                                      const Detector &d,
                                      const Scan &s) {
    for (std::size_t i = 0; i < obj.size(); ++i) {
      obj[i].update_centroid_mm(d, s);
    }
  }

  /** @returns The resolution of each observation */
  af::shared<double> observation_resolution(const af::const_ref<Observation> &obj,
                                            const BeamBase &b,
                                            const Detector &d) {
    af::shared<double> result(obj.size(), af::init_functor_null<double>());
    for (std::size_t i = 0; i < obj.size(); ++i) {
      result[i] = obj[i].resolution(b, d);
    }
    return result;
  }

  /**
   * A class to create a string from observations for pickling
   */
  struct observation_to_string : pickle_double_buffered::to_string {
    using pickle_double_buffered::to_string::operator<<;

    /** Initialise with the version */
    observation_to_string() {
      unsigned int version = 1;
      *this << version;
    }

    /** Convert a single observation to string */
    observation_to_string &operator<<(const Observation &val) {
      *this << val.panel << val.centroid.px.position[0] << val.centroid.px.position[1]
            << val.centroid.px.position[2] << val.centroid.px.variance[0]
            << val.centroid.px.variance[1] << val.centroid.px.variance[2]
            << val.centroid.px.std_err_sq[0] << val.centroid.px.std_err_sq[1]
            << val.centroid.px.std_err_sq[2] << val.centroid.mm.position[0]
            << val.centroid.mm.position[1] << val.centroid.mm.position[2]
            << val.centroid.mm.variance[0] << val.centroid.mm.variance[1]
            << val.centroid.mm.variance[2] << val.centroid.mm.std_err_sq[0]
            << val.centroid.mm.std_err_sq[1] << val.centroid.mm.std_err_sq[2]
            << val.intensity.observed.value << val.intensity.observed.variance
            << val.intensity.corrected.value << val.intensity.corrected.variance;

      return *this;
    }
  };

  /**
   * A class to convert strings to observations for unpicklings
   */
  struct observation_from_string : pickle_double_buffered::from_string {
    using pickle_double_buffered::from_string::operator>>;

    /** Initialise with the string */
    observation_from_string(const char *str_ptr)
        : pickle_double_buffered::from_string(str_ptr) {
      *this >> version;
      DIALS_ASSERT(version == 1);
    }

    /** Convert a string to a single observation instance */
    observation_from_string &operator>>(Observation &val) {
      *this >> val.panel >> val.centroid.px.position[0] >> val.centroid.px.position[1]
        >> val.centroid.px.position[2] >> val.centroid.px.variance[0]
        >> val.centroid.px.variance[1] >> val.centroid.px.variance[2]
        >> val.centroid.px.std_err_sq[0] >> val.centroid.px.std_err_sq[1]
        >> val.centroid.px.std_err_sq[2] >> val.centroid.mm.position[0]
        >> val.centroid.mm.position[1] >> val.centroid.mm.position[2]
        >> val.centroid.mm.variance[0] >> val.centroid.mm.variance[1]
        >> val.centroid.mm.variance[2] >> val.centroid.mm.std_err_sq[0]
        >> val.centroid.mm.std_err_sq[1] >> val.centroid.mm.std_err_sq[2]
        >> val.intensity.observed.value >> val.intensity.observed.variance
        >> val.intensity.corrected.value >> val.intensity.corrected.variance;

      return *this;
    }

    unsigned int version;
  };

  void export_flex_observation() {
    scitbx::af::boost_python::flex_wrapper<Observation, return_internal_reference<> >::
      plain("observation")
        .def("__init__", make_constructor(&init_from_panel))
        .def("__init__", make_constructor(&init_from_centroid))
        .def("__init__", make_constructor(&init_from_intensity))
        .def("__init__", make_constructor(&init_from_centroid_and_intensity))
        .def("__init__", make_constructor(&init_from_panel_and_centroid))
        .def("__init__", make_constructor(&init_from_panel_and_intensity))
        .def("__init__", make_constructor(&init_from_panel_centroid_and_intensity))
        .def("__init__", make_constructor(&init_from_single_panel_and_centroid))
        .def("__init__", make_constructor(&init_from_single_panel_and_intensity))
        .def("__init__",
             make_constructor(&init_from_single_panel_centroid_and_intensity))
        .def("panels", &observation_get_panels)
        .def("panels", &observation_set_panels)
        .def("centroids", &observation_get_centroids)
        .def("centroids", &observation_set_centroids)
        .def("intensities", &observation_get_intensities)
        .def("intensities", &observation_set_intensities)
        .def("update_centroid_mm",
             &observation_update_centroid_mm,
             (boost::python::arg("detector"), boost::python::arg("scan")))
        .def("resolution",
             &observation_resolution,
             (boost::python::arg("beam"), boost::python::arg("detector")))
        .def_pickle(flex_pickle_double_buffered<Observation,
                                                observation_to_string,
                                                observation_from_string>());
  }

}}}  // namespace dials::af::boost_python
