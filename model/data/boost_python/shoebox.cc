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
#include <dials/config.h>

namespace dials { namespace model { namespace boost_python {

  using namespace boost::python;
  using scitbx::af::flex;
  using scitbx::af::flex_int;

  /** Set the data array as a flex array */
  template <typename FloatType>
  void set_data(Shoebox<FloatType> &obj,
      typename flex<FloatType>::type &data) {
    DIALS_ASSERT(data.accessor().all().size() == 3);
    obj.data = af::versa<FloatType, af::c_grid<3> >(
      data.handle(), af::c_grid<3>(data.accessor()));
  }

  /** Set the mask array as a flex array */
  template <typename FloatType>
  void set_mask(Shoebox<FloatType> &obj, flex_int &mask) {
    DIALS_ASSERT(mask.accessor().all().size() == 3);
    obj.mask = af::versa<int, af::c_grid<3> >(
      mask.handle(), af::c_grid<3>(mask.accessor()));
  }

  /** Set the bgrd array as a flex array */
  template <typename FloatType>
  void set_background(Shoebox<FloatType> &obj,
      typename flex<FloatType>::type &background) {
    DIALS_ASSERT(background.accessor().all().size() == 3);
    obj.background = af::versa<FloatType, af::c_grid<3> >(
      background.handle(), af::c_grid<3>(background.accessor()));
  }

  template <typename FloatType>
  class_< Shoebox<FloatType> > shoebox_wrapper(const char *name)
  {
    typedef Shoebox<FloatType> shoebox_type;

    class_<shoebox_type> shoebox(name);
    shoebox
      .def(init<const shoebox_type&>())
      .def(init<std::size_t, const int6&>())
      .def(init<const int6&>())
      .add_property("data",
        make_getter(&shoebox_type::data,
          return_value_policy<return_by_value>()),
        &set_data<FloatType>)
      .add_property("mask",
        make_getter(&shoebox_type::mask,
          return_value_policy<return_by_value>()),
        &set_mask<FloatType>)
      .add_property("background",
        make_getter(&shoebox_type::background,
          return_value_policy<return_by_value>()),
        &set_background<FloatType>)
      .add_property("bbox",
        make_getter(&shoebox_type::bbox,
          return_value_policy<return_by_value>()),
        make_setter(&shoebox_type::bbox,
          return_value_policy<return_by_value>()))
      .def_readwrite("panel", &shoebox_type::panel)
      .def("allocate", &shoebox_type::allocate)
      .def("allocate", &shoebox_type::allocate_with_value)
      .def("deallocate", &shoebox_type::deallocate)
      .def("xoffset", &shoebox_type::xoffset)
      .def("yoffset", &shoebox_type::yoffset)
      .def("zoffset", &shoebox_type::zoffset)
      .def("offset", &shoebox_type::offset)
      .def("xsize", &shoebox_type::xsize)
      .def("ysize", &shoebox_type::ysize)
      .def("zsize", &shoebox_type::zsize)
      .def("size", &shoebox_type::size)
      .def("is_consistent", &shoebox_type::is_consistent)
      .def("is_bbox_within_image_volume",
        &shoebox_type::is_bbox_within_image_volume, (
          arg("image_size"), arg("scan_range")))
      .def("does_bbox_contain_bad_pixels",
        &shoebox_type::does_bbox_contain_bad_pixels, (
          arg("mask")))
      .def("count_mask_values",
        &shoebox_type::count_mask_values)
      .def("centroid_all",
        &shoebox_type::centroid_all)
      .def("centroid_masked",
        &shoebox_type::centroid_masked)
      .def("centroid_valid",
        &shoebox_type::centroid_valid)
      .def("centroid_foreground",
        &shoebox_type::centroid_foreground)
      .def("centroid_strong",
        &shoebox_type::centroid_strong)
      .def("centroid_all_minus_background",
        &shoebox_type::centroid_all_minus_background)
      .def("centroid_masked_minus_background",
        &shoebox_type::centroid_masked_minus_background)
      .def("centroid_valid_minus_background",
        &shoebox_type::centroid_valid_minus_background)
      .def("centroid_foreground_minus_background",
        &shoebox_type::centroid_foreground_minus_background)
      .def("centroid_strong_minus_background",
        &shoebox_type::centroid_strong_minus_background)
      .def("summed_intensity_all",
        &shoebox_type::summed_intensity_all)
      .def("summed_intensity_masked",
        &shoebox_type::summed_intensity_masked)
      .def("summed_intensity_valid",
        &shoebox_type::summed_intensity_valid)
      .def("summed_intensity_foreground",
        &shoebox_type::summed_intensity_foreground)
      .def("summed_intensity_strong",
        &shoebox_type::summed_intensity_strong)
      .def("__eq__", &shoebox_type::operator==)
      .def("__ne__", &shoebox_type::operator!=);

    return shoebox;
  }

  void export_shoebox()
  {
    shoebox_wrapper<ProfileFloatType>("Shoebox");
  }

}}} // namespace dials::model::boost_python
