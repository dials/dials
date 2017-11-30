/*
 * parallel_integrator_ext.cc
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
#include <dials/algorithms/integration/parallel_integrator.h>
#include <dials/algorithms/integration/algorithms.h>

using namespace boost::python;

namespace dials { namespace algorithms { namespace boost_python {

  /**
   * Get the reference profile data
   * @param data The reference profile structure
   * @index The index of the reference profile
   * @returns The data
   */
  af::versa< double, af::c_grid<3> > ReferenceProfileData_data(
      const ReferenceProfileData &data,
      std::size_t index) {
    af::versa< double, af::c_grid<3> > result(data.data(index).accessor());
    std::copy(
        data.data(index).begin(),
        data.data(index).end(),
        result.begin());
    return result;
  }

  /**
   * Get the reference profile mask
   * @param data The reference profile structure
   * @index The index of the reference profile
   * @returns The mask
   */
  af::versa< bool, af::c_grid<3> > ReferenceProfileData_mask(
      const ReferenceProfileData &data,
      std::size_t index) {
    af::versa< bool, af::c_grid<3> > result(data.mask(index).accessor());
    std::copy(
        data.mask(index).begin(),
        data.mask(index).end(),
        result.begin());
    return result;
  }

  /**
   * A class to pickle the reference profile data
   */
  struct ReferenceProfileDataPickleSuite : boost::python::pickle_suite {

    static
    boost::python::tuple getstate(const ReferenceProfileData& obj) {
      std::size_t version = 1;
      return boost::python::make_tuple(
          version,
          ReferenceProfileDataPickleSuite::get_data_list(obj),
          ReferenceProfileDataPickleSuite::get_mask_list(obj));
    }

    static
    void setstate(ReferenceProfileData& obj, boost::python::tuple state) {
      typedef af::const_ref<double, af::c_grid<3> > data_type;
      typedef af::const_ref<bool, af::c_grid<3> > mask_type;
      DIALS_ASSERT(boost::python::len(state) == 3);
      std::size_t version = boost::python::extract<std::size_t>(state[0])();
      DIALS_ASSERT(version == 1);
      DIALS_ASSERT(boost::python::len(state[1]) == boost::python::len(state[2]));
      for (std::size_t i = 0; i < boost::python::len(state[1]); ++i) {
        data_type data = boost::python::extract<data_type>(state[1][i])();
        mask_type mask = boost::python::extract<mask_type>(state[2][i])();
        obj.append(data, mask);
      }
    }

    static
    boost::python::list get_data_list(const ReferenceProfileData &obj) {
      boost::python::list result;
      for (std::size_t i = 0; i < obj.size(); ++i) {
        result.append(ReferenceProfileData_data(obj, i));
      }
      return result;
    }

    static
    boost::python::list get_mask_list(const ReferenceProfileData &obj) {
      boost::python::list result;
      for (std::size_t i = 0; i < obj.size(); ++i) {
        result.append(ReferenceProfileData_mask(obj, i));
      }
      return result;
    }

  };


  /**
   * A class to pickle the reference profile data
   */
  struct GaussianRSReferenceProfileDataPickleSuite : boost::python::pickle_suite {
    static
    boost::python::tuple getinitargs(const GaussianRSReferenceProfileData &obj) {
      return boost::python::make_tuple(
          obj.reference(),
          obj.sampler(),
          obj.spec());
    }
  };

  /**
   * A class to pickle the reference profile data
   */
  struct GaussianRSMultiCrystalReferenceProfileDataPickleSuite : boost::python::pickle_suite {

    static
    boost::python::tuple getstate(const GaussianRSMultiCrystalReferenceProfileData& obj) {
      std::size_t version = 1;
      return boost::python::make_tuple(
          version,
          GaussianRSMultiCrystalReferenceProfileDataPickleSuite::get_data_list(obj));
    }

    static
    void setstate(GaussianRSMultiCrystalReferenceProfileData& obj, boost::python::tuple state) {
      typedef GaussianRSReferenceProfileData data_type;
      DIALS_ASSERT(boost::python::len(state) == 2);
      std::size_t version = boost::python::extract<std::size_t>(state[0])();
      DIALS_ASSERT(version == 1);
      for (std::size_t i = 0; i < boost::python::len(state[1]); ++i) {
        data_type data = boost::python::extract<data_type>(state[1][i])();
        obj.append(data);
      }
    }

    static
    boost::python::list get_data_list(const GaussianRSMultiCrystalReferenceProfileData &obj) {
      boost::python::list result;
      for (std::size_t i = 0; i < obj.size(); ++i) {
        result.append(obj[i]);
      }
      return result;
    }
  };

  /**
   * Define the call operator of the IntensityCalculatorIface class
   * @param self The interface object
   * @param reflection The reflection to integrate
   * @param adjacent_reflections The reflections to integrator
   */
  void IntensityCalculatorIface_call(
        const IntensityCalculatorIface &self,
        af::Reflection &reflection,
        boost::python::object adjacent_reflections) {
    std::vector<af::Reflection> ar;
    for (std::size_t i = 0; i < boost::python::len(adjacent_reflections); ++i) {
      ar.push_back(boost::python::extract<af::Reflection>(adjacent_reflections[i])());
    }
    self(reflection, ar);
  }


  /**
   * Export the interfaces
   */
  void export_algorithm_interfaces() {

    class_<MaskCalculatorIface, boost::noncopyable>("MaskCalculatorIface", no_init)
      .def("__call__", &MaskCalculatorIface::operator())
      ;

    class_<BackgroundCalculatorIface, boost::noncopyable>("BackgroundCalculatorIface", no_init)
      .def("__call__", &BackgroundCalculatorIface::operator())
      ;

    class_<IntensityCalculatorIface, boost::noncopyable>("IntensityCalculatorIface", no_init)
      .def("__call__", &IntensityCalculatorIface_call)
      ;

  }

  /**
   * Export algorithms
   */
  void export_algorithms() {

    // Export GausianRSMaskCalculator
    class_<GaussianRSMaskCalculator,
           bases<MaskCalculatorIface> >("MaskCalculator", no_init)
      .def(init<const BeamBase&,
                const Detector&,
                const Goniometer&,
                const Scan&,
                double,
                double>())
      ;


    // Export GaussianRSMultiCrystalMaskCalculator
    class_<GaussianRSMultiCrystalMaskCalculator,
           bases<MaskCalculatorIface> >("MultiCrystalMaskCalculator")
      .def("append", &GaussianRSMultiCrystalMaskCalculator::append)
      ;


    // Export Simple background calculator
    class_<SimpleBackgroundCalculator,
           bases<BackgroundCalculatorIface> >("SimpleBackgroundCalculator", no_init)
      .def(init<boost::shared_ptr<background::Modeller>,
                boost::shared_ptr<background::OutlierRejector>,
                std::size_t>())
      ;

    // Export GLM background calculator
    class_<GLMBackgroundCalculator,
           bases<BackgroundCalculatorIface> >("GLMBackgroundCalculator", no_init)
      .def(init<GLMBackgroundCreator::Model,
                double,
                std::size_t,
                std::size_t>((
        arg("model"),
        arg("tuning_constant"),
        arg("max_iter"),
        arg("min_pixels"))))
      ;

    // Export GModel background calculator
    class_<GModelBackgroundCalculator,
           bases<BackgroundCalculatorIface> >("GModelBackgroundCalculator", no_init)
      .def(init<
                boost::shared_ptr<BackgroundModel>,
                bool,
                double,
                std::size_t,
                std::size_t>((
        arg("model"),
        arg("robust"),
        arg("tuning_constant"),
        arg("max_iter"),
        arg("min_pixels"))))
      ;


    // Export the reference data structure
    class_<ReferenceProfileData>("ReferenceProfileData")
      .def("append", &ReferenceProfileData::append)
      .def("data", &ReferenceProfileData_data)
      .def("mask", &ReferenceProfileData_mask)
      .def("__len__", &ReferenceProfileData::size)
      .def_pickle(ReferenceProfileDataPickleSuite())
      ;


    // Export GaussianRSReferenceProfileData
    class_<GaussianRSReferenceProfileData>("GaussianRSReferenceProfileData", no_init)
      .def(init<
          const ReferenceProfileData&,
          boost::shared_ptr<SamplerIface>,
          const TransformSpec&>())
      .def("reference",
          &GaussianRSReferenceProfileData::reference,
          return_internal_reference<>())
      .def("sampler",
          &GaussianRSReferenceProfileData::sampler)
      .def("spec",
          &GaussianRSReferenceProfileData::spec,
          return_internal_reference<>())
      .def_pickle(GaussianRSReferenceProfileDataPickleSuite())
      ;

    // Export GaussianRSMultiCrystalReferenceProfileData
    class_<GaussianRSMultiCrystalReferenceProfileData>("GaussianRSMultiCrystalReferenceProfileData")
      .def("append", &GaussianRSMultiCrystalReferenceProfileData::append)
      .def("__getitem__",
          &GaussianRSMultiCrystalReferenceProfileData::operator[],
          return_internal_reference<>())
      .def("__len__", &GaussianRSMultiCrystalReferenceProfileData::size)
      .def_pickle(GaussianRSMultiCrystalReferenceProfileDataPickleSuite())
      ;


    // Export NullIntensityCalculator
    class_<NullIntensityCalculator,
           bases<IntensityCalculatorIface> >("NullIntensityCalculator")
      ;


    // Export GaussianRSIntensityCalculator
    class_<GaussianRSIntensityCalculator,
           bases<IntensityCalculatorIface> >("GaussianRSIntensityCalculator", no_init)
      .def(init<
          const GaussianRSMultiCrystalReferenceProfileData &,
          bool,
          bool>())
      ;

  }

  /**
   * Export integrator
   */
  void export_integrator() {

    class_<ParallelIntegrator>("Integrator", no_init)
      .def(init<
          const af::reflection_table&,
          ImageSweep,
          const MaskCalculatorIface&,
          const BackgroundCalculatorIface&,
          const IntensityCalculatorIface&,
          std::size_t,
          bool>((
              arg("reflections"),
              arg("imageset"),
              arg("compute_mask"),
              arg("compute_background"),
              arg("compute_intensity"),
              arg("nthreads") = 1,
              arg("debug") = false)))
      .def("reflections",
          &ParallelIntegrator::reflections)
      ;

  }

  BOOST_PYTHON_MODULE(dials_algorithms_integration_parallel_integrator_ext)
  {
    export_algorithm_interfaces();
    export_algorithms();
    export_integrator();
  }

}}} // namespace dials_scratch::examples::boost_python
