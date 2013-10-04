/*
 * extractor.cc
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
#include <dials/algorithms/shoebox/extractor.h>

namespace dials { namespace algorithms { namespace shoebox {
  namespace boost_python {

  using namespace boost::python;

  /**
   * A constructor for multiple panels
   */
  static
  Extractor* make_for_multi_panel(
      const af::const_ref<std::size_t> &panels,
      const af::const_ref<int6> &bboxes,
      const boost::python::tuple &mask,
      const boost::python::tuple &gain,
      const boost::python::tuple &dark) {

    typedef af::versa< bool, af::flex_grid<> > versa_bool2;
    typedef af::versa< double, af::flex_grid<> > versa_double2;

    // Check the number of panels is consistent
    std::size_t npanels = boost::python::len(mask);
    DIALS_ASSERT(npanels > 0);
    DIALS_ASSERT(boost::python::len(gain) == npanels);
    DIALS_ASSERT(boost::python::len(dark) == npanels);    

    // Get the shape of all the panels into an array and check the size of
    // each is consistent. Get the total array size needed
    af::shared<int2> shape(npanels, af::init_functor_null<int2>());
    std::size_t size = 0;
    for (std::size_t i = 0; i < npanels; ++i) {
      versa_bool2 m = boost::python::extract<versa_bool2>(mask[i]);
      versa_double2 g = boost::python::extract<versa_double2>(gain[i]);
      versa_double2 d = boost::python::extract<versa_double2>(dark[i]);
      DIALS_ASSERT(m.accessor().all().size() == 2);
      DIALS_ASSERT(m.accessor().all().all_gt(0));
      DIALS_ASSERT(m.accessor().all().all_eq(g.accessor().all()));
      DIALS_ASSERT(m.accessor().all().all_eq(d.accessor().all()));
      shape[i][0] = m.accessor().all()[0];
      shape[i][1] = m.accessor().all()[1];
      size += m.size();
    }
    
    // Construct the 1D arrays
    af::shared<bool> mask_1d(size, af::init_functor_null<bool>());
    af::shared<double> gain_1d(size, af::init_functor_null<double>());
    af::shared<double> dark_1d(size, af::init_functor_null<double>());
    std::size_t k = 0;
    for (std::size_t i = 0; i < npanels; ++i) {
      versa_bool2 m = boost::python::extract<versa_bool2>(mask[i]);
      versa_double2 g = boost::python::extract<versa_double2>(gain[i]);
      versa_double2 d = boost::python::extract<versa_double2>(dark[i]);
      for (std::size_t j = 0; j < m.size(); ++j, ++k) {
        mask_1d[k] = m[j];
        gain_1d[k] = g[j];
        dark_1d[k] = d[j];
      }  
    }
    
    // Return the new extractor object
    return new Extractor(panels, bboxes, mask_1d, gain_1d, dark_1d, shape);
  }

  /**
   * A constructor for single panels
   */
  static
  Extractor* make_for_single_panel(
      const af::const_ref<int6> bboxes, 
      af::versa< bool, scitbx::af::flex_grid<> > mask,
      af::versa< double, scitbx::af::flex_grid<> > gain,
      af::versa< double, scitbx::af::flex_grid<> > dark) {
    af::c_grid<2> accessor(mask.accessor());
    return new Extractor(
      bboxes, 
      af::versa< bool, af::c_grid<2> >(mask.handle(), accessor),
      af::versa< double, af::c_grid<2> >(gain.handle(), accessor),
      af::versa< double, af::c_grid<2> >(dark.handle(), accessor));
  }

  void export_extractor()
  {
    class_ <Extractor> ("Extractor", no_init)
      .def(init<const af::const_ref<std::size_t>&,
                const af::const_ref<int6>&,
                const af::shared<bool>&,
                const af::shared<double>&,
                const af::shared<double>&,
                const af::shared<int2>&>((
          arg("panels"),
          arg("bboxes"),
          arg("mask"),
          arg("gain"),
          arg("dark"),
          arg("shape"))))
      .def("__init__", make_constructor(
        &make_for_multi_panel, 
        default_call_policies(), (
          arg("panels"),
          arg("bboxes"),
          arg("mask"),
          arg("gain"),
          arg("dark"))))
      .def("__init__", make_constructor(
        &make_for_single_panel, 
        default_call_policies(), (
          arg("bboxes"),
          arg("mask"),
          arg("gain"),
          arg("dark"))))
      .def("add_image", &Extractor::add_image, (
        arg("panel"), 
        arg("frame"), 
        arg("image")))
      .def("indices", &Extractor::indices)
      .def("shoeboxes", &Extractor::shoeboxes);
  }

}}}} // namespace = dials::algorithms::shoebox::boost_python
