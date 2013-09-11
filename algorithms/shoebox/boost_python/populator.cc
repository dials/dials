/*
 * populator.cc
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
#include <dials/algorithms/shoebox/populator.h>

namespace dials { namespace algorithms { namespace shoebox {
  namespace boost_python {

  using namespace boost::python;

  Populator* make_populator(
      af::shared<Reflection> reflections, 
      af::versa< bool, scitbx::af::flex_grid<> > mask,
      af::versa< double, scitbx::af::flex_grid<> > gain,
      af::versa< double, scitbx::af::flex_grid<> > dark) {
    DIALS_ASSERT(mask.accessor().all().all_eq(gain.accessor().all()));
    DIALS_ASSERT(mask.accessor().all().all_eq(dark.accessor().all()));
    af::c_grid<2> accessor(mask.accessor());
    return new Populator(
      reflections, 
      af::versa< bool, af::c_grid<2> >(mask.handle(), accessor),
      af::versa< double, af::c_grid<2> >(gain.handle(), accessor),
      af::versa< double, af::c_grid<2> >(dark.handle(), accessor));
  }

  void export_populator()
  {
    class_ <Populator> ("Populator", no_init)
      .def("__init__", make_constructor(
        &make_populator, 
        default_call_policies(), (
          arg("reflection_list"),
          arg("mask"),
          arg("gain_map"),
          arg("dark_map"))))
      .def("add_image", &Populator::add_image, (
        arg("image"), arg("index")))
      .def("image_mask", &Populator::image_mask, (
        arg("index"), arg("kernal_size")))
      .def("indices", &Populator::indices);
  }

}}}} // namespace = dials::algorithms::shoebox::boost_python
