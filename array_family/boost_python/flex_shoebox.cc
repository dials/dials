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
  using scitbx::af::int2;
  using scitbx::af::int6;
  using scitbx::af::small;
  using scitbx::vec3;
  using dials::model::Shoebox;
  using dials::model::Valid;
  using dials::model::Foreground;
  using dials::algorithms::LabelImageStack;

  /**
   * Construct an array of shoebxoes from a spot labelling class
   */
  scitbx::af::flex<Shoebox>::type* from_labels(const LabelImageStack &label) {

    // Get the stuff from the label struct
    scitbx::af::shared<int> labels = label.labels();
    scitbx::af::shared<int> values = label.values();
    scitbx::af::shared< vec3<int> > coords = label.coords();

    // Get the number of labels and allocate the array
    std::size_t num = scitbx::af::max(labels.const_ref()) + 1;
    scitbx::af::shared<Shoebox> result(num);
    
    // Initialise the bboxes
    int xsize = label.size()[1];
    int ysize = label.size()[0];
    int zsize = label.num_images();
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i].bbox[0] = xsize; result[i].bbox[1] = 0;
      result[i].bbox[2] = ysize; result[i].bbox[3] = 0;
      result[i].bbox[4] = zsize; result[i].bbox[5] = 0;
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
      result[l].mask(kk,jj,ii) = Valid | Foreground;
    }  

    // Return the array
    return new scitbx::af::flex<Shoebox>::type(
      result, scitbx::af::flex_grid<>(num));
  }  
  
  /**
   * Check if the arrays are consistent
   */
  scitbx::af::flex_bool is_consistent(
      const scitbx::af::const_ref<Shoebox> &a) {
    scitbx::af::flex_bool result(a.size());
    for (std::size_t i = 0; i < a.size(); ++i) {
      result[i] = a[i].is_consistent();
    }
    return result;
  }
  
  /**
   * Check if the bounding box has points outside the image range.
   */
  scitbx::af::flex_bool is_bbox_within_image_volume(
      const scitbx::af::const_ref<Shoebox> &a,
      small<long,10> image_size, int2 scan_range) {
    scitbx::af::flex_bool result(a.size());
    for (std::size_t i = 0; i < a.size(); ++i) {
      result[i] = a[i].is_bbox_within_image_volume(image_size, scan_range);
    }
    return result;
  }
  
  /**
   * Check if the bounding box has points that cover bad pixels
   */
  scitbx::af::flex_bool does_bbox_contain_bad_pixels(
      const scitbx::af::const_ref<Shoebox> &a, 
      const scitbx::af::flex_bool &mask) {
    scitbx::af::flex_bool result(a.size());
    for (std::size_t i = 0; i < a.size(); ++i) {
      result[i] = a[i].does_bbox_contain_bad_pixels(mask);
    }
    return result;
  }
  
  /**
   * Count the number of mask pixels with the given code
   */
  scitbx::af::flex_int count_mask_values(
      const scitbx::af::const_ref<Shoebox> &a,
      int code) {
    scitbx::af::flex_int result(a.size());
    for (std::size_t i = 0; i < a.size(); ++i) {
      result[i] = a[i].count_mask_values(code);
    }
    return result;
  }

  /**
   * Get the bounding boxes
   */
  scitbx::af::flex<int6>::type bounding_boxes(
      const scitbx::af::const_ref<Shoebox> &a) {
    scitbx::af::flex<int6>::type result(a.size());
    for (std::size_t i = 0; i < a.size(); ++i) {
      result[i] = a[i].bbox;
    }
    return result;
  }

  void export_flex_shoebox()
  {
    scitbx::af::boost_python::flex_wrapper <
        Shoebox, return_internal_reference<> >::plain("shoebox")
      .def("__init__", make_constructor(from_labels))
      .def("is_consistent", &is_consistent)
      .def("bounding_boxes", &bounding_boxes)
      .def("count_mask_values", &count_mask_values)
      .def("is_bbox_within_image_volume", &is_bbox_within_image_volume, (
        arg("image_size"), 
        arg("scan_range")))
      .def("does_bbox_contain_bad_pixels", &does_bbox_contain_bad_pixels, (
        arg("mask")));
  }

}}} // namespace dials::af::boost_python
