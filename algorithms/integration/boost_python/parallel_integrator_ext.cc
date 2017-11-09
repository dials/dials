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
                std::size_t>())
      ;

    // Export GModel background calculator
    class_<GModelBackgroundCalculator,
           bases<BackgroundCalculatorIface> >("GModelBackgroundCalculator", no_init)
      .def(init<
                boost::shared_ptr<BackgroundModel>,
                bool,
                double,
                std::size_t,
                std::size_t>())
      ;


    // Export the reference data structure
    class_<Reference>("Reference")
      .def("append", &Reference::append)
      ;


    // Export GaussianRSIntensityCalculatorData
    class_<GaussianRSIntensityCalculatorData>("GaussianRSIntensityCalculatorData", no_init)
      .def(init<
          const Reference&,
          const CircleSampler&,
          const TransformSpec&>())
      ;

    // Export GaussianRSMultiCrystalIntensityCalculatorData
    class_<GaussianRSMultiCrystalIntensityCalculatorData>("GaussianRSMultiCrystalIntensityCalculatorData")
      .def("append", &GaussianRSMultiCrystalIntensityCalculatorData::append)
      ;


    // Export NullIntensityCalculator
    class_<NullIntensityCalculator,
           bases<IntensityCalculatorIface> >("NullIntensityCalculator")
      ;


    // Export GaussianRSIntensityCalculator
    class_<GaussianRSIntensityCalculator,
           bases<IntensityCalculatorIface> >("GaussianRSIntensityCalculator", no_init)
      .def(init<
          const GaussianRSMultiCrystalIntensityCalculatorData &,
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
