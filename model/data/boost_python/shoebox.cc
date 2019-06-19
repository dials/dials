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
  void set_data(Shoebox<FloatType> &obj, typename flex<FloatType>::type &data) {
    DIALS_ASSERT(data.accessor().all().size() == 3);
    obj.data = af::versa<FloatType, af::c_grid<3> >(data.handle(),
                                                    af::c_grid<3>(data.accessor()));
  }

  /** Set the mask array as a flex array */
  template <typename FloatType>
  void set_mask(Shoebox<FloatType> &obj, flex_int &mask) {
    DIALS_ASSERT(mask.accessor().all().size() == 3);
    obj.mask =
      af::versa<int, af::c_grid<3> >(mask.handle(), af::c_grid<3>(mask.accessor()));
  }

  /** Set the bgrd array as a flex array */
  template <typename FloatType>
  void set_background(Shoebox<FloatType> &obj,
                      typename flex<FloatType>::type &background) {
    DIALS_ASSERT(background.accessor().all().size() == 3);
    obj.background = af::versa<FloatType, af::c_grid<3> >(
      background.handle(), af::c_grid<3>(background.accessor()));
  }

  /** Get a list of shoebox pixel coordinates. */
  template <typename ShoeboxType>
  af::shared<vec3<double> > coords(const ShoeboxType &self) {
    af::shared<vec3<double> > result;
    for (std::size_t k = 0; k < self.zsize(); ++k) {
      for (std::size_t j = 0; j < self.ysize(); ++j) {
        for (std::size_t i = 0; i < self.xsize(); ++i) {
          double x = i + self.xoffset() + 0.5;
          double y = j + self.yoffset() + 0.5;
          double z = k + self.zoffset() + 0.5;
          result.push_back(vec3<double>(x, y, z));
        }
      }
    }
    return result;
  }

  /** Get a list of shoebox pixel coordinates. */
  template <typename ShoeboxType>
  af::shared<vec3<double> > coords_with_mask(
    const ShoeboxType &self,
    const af::const_ref<bool, af::c_grid<3> > mask) {
    DIALS_ASSERT(mask.accessor().all_eq(self.size()));
    af::shared<vec3<double> > result;
    for (std::size_t k = 0; k < self.zsize(); ++k) {
      for (std::size_t j = 0; j < self.ysize(); ++j) {
        for (std::size_t i = 0; i < self.xsize(); ++i) {
          if (mask(k, j, i)) {
            double x = i + self.xoffset() + 0.5;
            double y = j + self.yoffset() + 0.5;
            double z = k + self.zoffset() + 0.5;
            result.push_back(vec3<double>(x, y, z));
          }
        }
      }
    }
    return result;
  }

  /** Get a list of shoebox pixel values. */
  template <typename ShoeboxType>
  af::shared<double> values(const ShoeboxType &self) {
    DIALS_ASSERT(self.is_consistent());
    af::shared<double> result;
    for (std::size_t k = 0; k < self.zsize(); ++k) {
      for (std::size_t j = 0; j < self.ysize(); ++j) {
        for (std::size_t i = 0; i < self.xsize(); ++i) {
          result.push_back(self.data(k, j, i));
        }
      }
    }
    return result;
  }

  /** Get a list of shoebox pixel values. */
  template <typename ShoeboxType>
  af::shared<double> values_with_mask(const ShoeboxType &self,
                                      const af::const_ref<bool, af::c_grid<3> > mask) {
    DIALS_ASSERT(mask.accessor().all_eq(self.size()));
    DIALS_ASSERT(self.is_consistent());
    af::shared<double> result;
    for (std::size_t k = 0; k < self.zsize(); ++k) {
      for (std::size_t j = 0; j < self.ysize(); ++j) {
        for (std::size_t i = 0; i < self.xsize(); ++i) {
          if (mask(k, j, i)) {
            result.push_back(self.data(k, j, i));
          }
        }
      }
    }
    return result;
  }

  /** Get a list of shoebox pixel coordinates. */
  template <typename ShoeboxType>
  af::shared<vec3<double> > beam_vectors(const ShoeboxType &self,
                                         const Detector &detector) {
    af::shared<vec3<double> > result;
    for (std::size_t k = 0; k < self.zsize(); ++k) {
      for (std::size_t j = 0; j < self.ysize(); ++j) {
        for (std::size_t i = 0; i < self.xsize(); ++i) {
          double x = i + self.xoffset() + 0.5;
          double y = j + self.yoffset() + 0.5;
          vec2<double> c(x, y);
          result.push_back(detector[self.panel].get_pixel_lab_coord(c));
        }
      }
    }
    return result;
  }

  /** Get a list of shoebox pixel coordinates. */
  template <typename ShoeboxType>
  af::shared<vec3<double> > beam_vectors_with_mask(
    const ShoeboxType &self,
    const Detector &detector,
    const af::const_ref<bool, af::c_grid<3> > mask) {
    DIALS_ASSERT(mask.accessor().all_eq(self.size()));
    af::shared<vec3<double> > result;
    for (std::size_t k = 0; k < self.zsize(); ++k) {
      for (std::size_t j = 0; j < self.ysize(); ++j) {
        for (std::size_t i = 0; i < self.xsize(); ++i) {
          if (mask(k, j, i)) {
            double x = i + self.xoffset() + 0.5;
            double y = j + self.yoffset() + 0.5;
            vec2<double> c(x, y);
            result.push_back(detector[self.panel].get_pixel_lab_coord(c));
          }
        }
      }
    }
    return result;
  }

  template <typename FloatType>
  class_<Shoebox<FloatType> > shoebox_wrapper(const char *name) {
    typedef Shoebox<FloatType> shoebox_type;

    class_<shoebox_type> shoebox(name);
    shoebox.def(init<const shoebox_type &>())
      .def(init<std::size_t, const int6 &>())
      .def(init<const int6 &>())
      .add_property(
        "data",
        make_getter(&shoebox_type::data, return_value_policy<return_by_value>()),
        &set_data<FloatType>)
      .add_property(
        "mask",
        make_getter(&shoebox_type::mask, return_value_policy<return_by_value>()),
        &set_mask<FloatType>)
      .add_property(
        "background",
        make_getter(&shoebox_type::background, return_value_policy<return_by_value>()),
        &set_background<FloatType>)
      .add_property(
        "bbox",
        make_getter(&shoebox_type::bbox, return_value_policy<return_by_value>()),
        make_setter(&shoebox_type::bbox, return_value_policy<return_by_value>()))
      .def_readwrite("panel", &shoebox_type::panel)
      .def_readwrite("flat", &shoebox_type::flat)
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
      .def("size_flat", &shoebox_type::size_flat)
      .def("is_consistent", &shoebox_type::is_consistent)
      .def("is_allocated", &shoebox_type::is_allocated)
      .def("is_bbox_within_image_volume",
           &shoebox_type::is_bbox_within_image_volume,
           (arg("image_size"), arg("scan_range")))
      .def("does_bbox_contain_bad_pixels",
           &shoebox_type::does_bbox_contain_bad_pixels,
           (arg("mask")))
      .def("count_mask_values", &shoebox_type::count_mask_values)
      .def("all_foreground_valid", &shoebox_type::all_foreground_valid)
      .def("centroid_all", &shoebox_type::centroid_all)
      .def("centroid_masked", &shoebox_type::centroid_masked)
      .def("centroid_valid", &shoebox_type::centroid_valid)
      .def("centroid_foreground", &shoebox_type::centroid_foreground)
      .def("centroid_strong", &shoebox_type::centroid_strong)
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
      .def("bayesian_intensity", &shoebox_type::bayesian_intensity)
      .def("summed_intensity", &shoebox_type::summed_intensity)
      .def("flatten", &shoebox_type::flatten)
      .def("coords", &coords<shoebox_type>)
      .def("coords", &coords_with_mask<shoebox_type>)
      .def("values", &values<shoebox_type>)
      .def("values", &values_with_mask<shoebox_type>)
      .def("beam_vectors", &beam_vectors<shoebox_type>)
      .def("beam_vectors", &beam_vectors_with_mask<shoebox_type>)
      .def("__eq__", &shoebox_type::operator==)
      .def("__ne__", &shoebox_type::operator!=);

    return shoebox;
  }

  void export_shoebox() {
    shoebox_wrapper<ProfileFloatType>("Shoebox");
  }

}}}  // namespace dials::model::boost_python
