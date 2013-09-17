/*
 * shoebox.cc
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
#include <scitbx/array_family/flex_types.h>
#include <dials/model/data/shoebox.h>

namespace dials { namespace model { namespace boost_python {

  using namespace boost::python;
  using scitbx::af::flex_double;
  using scitbx::af::flex_int;
   
  /** Set the data array as a flex array */
  static
  void set_data(Shoebox &obj, flex_double &data) {
    DIALS_ASSERT(data.accessor().all().size() == 3);
    obj.data = af::versa<double, af::c_grid<3> >(
      data.handle(), af::c_grid<3>(data.accessor()));
  }

  /** Set the mask array as a flex array */
  static
  void set_mask(Shoebox &obj, flex_int &mask) {
    DIALS_ASSERT(mask.accessor().all().size() == 3);
    obj.mask = af::versa<int, af::c_grid<3> >(
      mask.handle(), af::c_grid<3>(mask.accessor()));
  }

  /** Set the bgrd array as a flex array */
  static
  void set_background(Shoebox &obj, flex_double &background) {
    DIALS_ASSERT(background.accessor().all().size() == 3);
    obj.background = af::versa<double, af::c_grid<3> >(
      background.handle(), af::c_grid<3>(background.accessor()));
  }

  void export_shoebox()
  {
    class_<Shoebox>("Shoebox")
      .def(init<const Shoebox&>())
      .def(init<std::size_t, const int6&>())
      .def(init<const int6&>())
      .add_property("data", 
        make_getter(&Shoebox::data, return_value_policy<return_by_value>()),
        &set_data)        
      .add_property("mask", 
        make_getter(&Shoebox::mask, return_value_policy<return_by_value>()),
        &set_mask)  
      .add_property("background", 
        make_getter(&Shoebox::background, 
          return_value_policy<return_by_value>()),
        &set_background)        
      .add_property("bbox", 
        make_getter(&Shoebox::bbox, return_value_policy<return_by_value>()),
        make_setter(&Shoebox::bbox, return_value_policy<return_by_value>())) 
      .def_readwrite("panel", &Shoebox::panel)       
      .def("allocate", &Shoebox::allocate)
      .def("deallocate", &Shoebox::deallocate)
      .def("xoffset", &Shoebox::xoffset)
      .def("yoffset", &Shoebox::yoffset)
      .def("zoffset", &Shoebox::zoffset)
      .def("offset", &Shoebox::offset)
      .def("xsize", &Shoebox::xsize)
      .def("ysize", &Shoebox::ysize)
      .def("zsize", &Shoebox::zsize)
      .def("size", &Shoebox::size)
      .def("is_consistent", &Shoebox::is_consistent)
      .def("is_bbox_within_image_volume", 
        &Shoebox::is_bbox_within_image_volume, (
          arg("image_size"), arg("scan_range")))
      .def("does_bbox_contain_bad_pixels",
        &Shoebox::does_bbox_contain_bad_pixels, (
          arg("mask")))
      .def("count_mask_values",
        &Shoebox::count_mask_values)
      .def("centroid_all",
        &Shoebox::centroid_all)
      .def("centroid_masked", 
        &Shoebox::centroid_masked)
      .def("centroid_valid", 
        &Shoebox::centroid_valid)
      .def("centroid_foreground", 
        &Shoebox::centroid_foreground)
      .def("centroid_strong", 
        &Shoebox::centroid_strong)
      .def("centroid_all_minus_background", 
        &Shoebox::centroid_all_minus_background)
      .def("centroid_masked_minus_background", 
        &Shoebox::centroid_masked_minus_background)
      .def("centroid_valid_minus_background", 
        &Shoebox::centroid_valid_minus_background)
      .def("centroid_foreground_minus_background", 
        &Shoebox::centroid_foreground_minus_background)
      .def("centroid_strong_minus_background", 
        &Shoebox::centroid_strong_minus_background)
      .def("summed_intensity_all",
        &Shoebox::summed_intensity_all)
      .def("summed_intensity_masked",
        &Shoebox::summed_intensity_masked)
      .def("summed_intensity_valid",
        &Shoebox::summed_intensity_valid)
      .def("summed_intensity_foreground",
        &Shoebox::summed_intensity_foreground)
      .def("summed_intensity_strong",
        &Shoebox::summed_intensity_strong)
      .def("__eq__", &Shoebox::operator==)
      .def("__ne__", &Shoebox::operator!=);
  }

}}} // namespace dials::model::boost_python
