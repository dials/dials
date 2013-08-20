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
#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <scitbx/array_family/ref_reductions.h>
#include <dials/model/data/shoebox.h>
#include <dials/algorithms/image/connected_components/connected_components.h>

namespace dials { namespace af { namespace boost_python {

  using namespace boost::python;
  using scitbx::vec3;
  using dials::model::Shoebox;
  using dials::algorithms::LabelImageStack;

  scitbx::af::flex<Shoebox>::type* from_labels(const LabelImageStack &label) {

    // Get the stuff from the label struct
    scitbx::af::shared<int> labels = label.labels();
    scitbx::af::shared<int> values = label.values();
    scitbx::af::shared< vec3<int> > coords = label.coords();

    // Get the number of labels and allocate the array
    std::size_t num = scitbx::af::max(labels.const_ref()) + 1;
    scitbx::af::shared<Shoebox> result(num);
    
    // Initialise the bboxes
    const int HIGH = 1000000;
    const int LOW = -1000000;
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i].bbox[0] = HIGH; result[i].bbox[1] = LOW;
      result[i].bbox[2] = HIGH; result[i].bbox[3] = LOW;
      result[i].bbox[4] = HIGH; result[i].bbox[5] = LOW;
    }

    // Set the shoeboxes
    for (std::size_t i = 0; i < labels.size(); ++i) {
      int l = labels[i];
      vec3<int> c = coords[i];
      if (c[2] <  result[l].bbox[0]) result[l].bbox[0] = c[2];
      if (c[2] >= result[l].bbox[1]) result[l].bbox[1] = c[2] + 1;
      if (c[1] <  result[l].bbox[2]) result[l].bbox[2] = c[1];
      if (c[1] >= result[l].bbox[3]) result[l].bbox[3] = c[1] + 1;
      if (c[0] <  result[l].bbox[4]) result[l].bbox[4] = c[0];
      if (c[0] >= result[l].bbox[5]) result[l].bbox[5] = c[0] + 1;
    }

    // Allocate all the arrays
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i].allocate();
    } 

    // Set all the mask and data points
    for (std::size_t i = 0; i < labels.size(); ++i) {
      int l = labels[i];
      double v = values[i];
      vec3<int> c = coords[i];
      int ii = c[2] - result[l].bbox[0];
      int jj = c[1] - result[l].bbox[2];
      int kk = c[0] - result[l].bbox[4];
      DIALS_ASSERT(ii >= 0 && jj >= 0 && kk >= 0);
      DIALS_ASSERT(ii < result[l].xsize());
      DIALS_ASSERT(jj < result[l].ysize());
      DIALS_ASSERT(kk < result[l].zsize());     
      result[l].data(kk,jj,ii) = (double)v;
      result[l].mask(kk,jj,ii) = 1;
    }  

    // Return the array
    return new scitbx::af::flex<Shoebox>::type(
      result, scitbx::af::flex_grid<>(num));
  }  

  void export_flex_shoebox()
  {
    scitbx::af::boost_python::flex_wrapper <Shoebox>::plain("shoebox")
      .def("__init__", make_constructor(from_labels));
  }

}}} // namespace dials::af::boost_python
