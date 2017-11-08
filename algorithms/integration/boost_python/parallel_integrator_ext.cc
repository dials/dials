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


  void export_algorithm_interfaces() {
    
    class_<MaskCalculatorIface, boost::noncopyable>("MaskCalculatorIface", no_init)
      ;

    class_<BackgroundCalculatorIface, boost::noncopyable>("BackgroundCalculatorIface", no_init)
      ;

    class_<IntensityCalculatorIface, boost::noncopyable>("IntensityCalculatorIface", no_init)
      ;

  }

  void export_algorithms() {

    class_<GaussianRSMaskCalculator, bases<MaskCalculatorIface> >("GaussianRSMaskCalculator", no_init)
      .def(init<const BeamBase&,
                const Detector&,
                const Goniometer&,
                const Scan&,
                double,
                double>())
      ;

    class_<GLMBackgroundCalculator, bases<BackgroundCalculatorIface> >("BackgroundCalculator", no_init)
      .def(init<GLMBackgroundCreator::Model,
                double,
                std::size_t,
                std::size_t>())
      ;


    class_<Reference>("Reference")
      .def("append", &Reference::append)
      ;

    class_<IntensityCalculator, bases<IntensityCalculatorIface> >("IntensityCalculator", no_init)
      .def(init<
          const Reference&,
          const CircleSampler&,
          const TransformSpec&,
          bool,
          bool>())
      ;

  }

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
