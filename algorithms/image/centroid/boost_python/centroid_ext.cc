/*
 * centroid_ext.cc
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
#include <dials/algorithms/image/centroid/centroid_points.h>
#include <dials/algorithms/image/centroid/centroid_image.h>
#include <dials/algorithms/image/centroid/centroid_masked_image.h>
#include <dials/algorithms/image/centroid/bias.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  template <typename FloatType, typename CoordType>
  class_<CentroidPoints<FloatType, CoordType> > centroid_points_wrapper(
    const char *name) {
    typedef CentroidPoints<FloatType, CoordType> CentroidPointsType;

    return class_<CentroidPointsType>(name, no_init)
      .def(init<const af::const_ref<FloatType> &, const af::const_ref<CoordType> &>(
        (arg("pixels"), arg("coord"))))
      .def("sum_pixels", &CentroidPointsType::sum_pixels)
      .def("sum_pixels_sq", &CentroidPointsType::sum_pixels_sq)
      .def("sum_pixels_coords", &CentroidPointsType::sum_pixels_coords)
      .def("sum_pixels_delta_sq", &CentroidPointsType::sum_pixels_delta_sq)
      .def("sum_pixels_delta_cross", &CentroidPointsType::sum_pixels_delta_cross)
      .def("mean", &CentroidPointsType::mean)
      .def("variance", &CentroidPointsType::variance)
      .def("unbiased_variance", &CentroidPointsType::unbiased_variance)
      .def("standard_error_sq", &CentroidPointsType::standard_error_sq)
      .def("unbiased_standard_error_sq",
           &CentroidPointsType::unbiased_standard_error_sq)
      .def("average_bias_estimate", &CentroidPointsType::average_bias_estimate)
      .def("mean_sq_error", &CentroidPointsType::mean_sq_error)
      .def("covariance_matrix", &CentroidPointsType::covariance_matrix);
  }

  template <typename FloatType>
  void centroid_image_2d_wrapper(const char *name) {
    typedef CentroidPoints<FloatType, vec2<double> > CentroidPointsType;

    class_<CentroidImage2d<FloatType>, bases<CentroidPointsType> >(name, no_init)
      .def(init<const af::const_ref<FloatType, af::c_grid<2> > &>((arg("image"))));
  }

  template <typename FloatType>
  void centroid_image_3d_wrapper(const char *name) {
    typedef CentroidPoints<FloatType, vec3<double> > CentroidPointsType;

    class_<CentroidImage3d<FloatType>, bases<CentroidPointsType> >(name, no_init)
      .def(init<const af::const_ref<FloatType, af::c_grid<3> > &>((arg("image"))));
  }

  template <typename FloatType>
  void centroid_masked_image_2d_wrapper(const char *name) {
    typedef CentroidPoints<FloatType, vec2<double> > CentroidPointsType;

    class_<CentroidMaskedImage2d<FloatType>, bases<CentroidPointsType> >(name, no_init)
      .def(
        init<const af::const_ref<FloatType, af::c_grid<2> > &,
             const af::const_ref<bool, af::c_grid<2> > &>((arg("image"), arg("mask"))));
  }

  template <typename FloatType>
  void centroid_masked_image_3d_wrapper(const char *name) {
    typedef CentroidPoints<FloatType, vec3<double> > CentroidPointsType;

    class_<CentroidMaskedImage3d<FloatType>, bases<CentroidPointsType> >(name, no_init)
      .def(
        init<const af::const_ref<FloatType, af::c_grid<3> > &,
             const af::const_ref<bool, af::c_grid<3> > &>((arg("image"), arg("mask"))));
  }

  template <typename FloatType, typename CoordType>
  CentroidPoints<FloatType, CoordType> centroid_points(
    const af::const_ref<FloatType> &pixels,
    const af::const_ref<CoordType> &coords) {
    return CentroidPoints<FloatType, CoordType>(pixels, coords);
  }

  template <typename FloatType>
  CentroidImage2d<FloatType> centroid_image_2d(
    const af::const_ref<FloatType, af::c_grid<2> > &image) {
    return CentroidImage2d<FloatType>(image);
  }

  template <typename FloatType>
  CentroidImage3d<FloatType> centroid_image_3d(
    const af::const_ref<FloatType, af::c_grid<3> > &image) {
    return CentroidImage3d<FloatType>(image);
  }

  template <typename FloatType>
  CentroidMaskedImage2d<FloatType> centroid_masked_image_2d(
    const af::const_ref<FloatType, af::c_grid<2> > &image,
    const af::const_ref<bool, af::c_grid<2> > &mask) {
    return CentroidMaskedImage2d<FloatType>(image, mask);
  }

  template <typename FloatType>
  CentroidMaskedImage3d<FloatType> centroid_masked_image_3d(
    const af::const_ref<FloatType, af::c_grid<3> > &image,
    const af::const_ref<bool, af::c_grid<3> > &mask) {
    return CentroidMaskedImage3d<FloatType>(image, mask);
  }

  template <typename FloatType>
  void centroid_suite() {
    def("centroid_points",
        &centroid_points<FloatType, vec2<double> >,
        (arg("pixels"), arg("coords")));

    def("centroid_points",
        &centroid_points<FloatType, vec3<double> >,
        (arg("pixels"), arg("coords")));

    def("centroid_image", &centroid_image_2d<FloatType>, (arg("image")));

    def("centroid_image", &centroid_image_3d<FloatType>, (arg("image")));

    def("centroid_image",
        &centroid_masked_image_2d<FloatType>,
        (arg("image"), arg("mask")));

    def("centroid_image",
        &centroid_masked_image_3d<FloatType>,
        (arg("image"), arg("mask")));
  }

  BOOST_PYTHON_MODULE(dials_algorithms_image_centroid_ext) {
    centroid_points_wrapper<float, vec2<double> >("CentroidPoints2dFloat");
    centroid_points_wrapper<float, vec3<double> >("CentroidPoints3dFloat");
    centroid_points_wrapper<double, vec2<double> >("CentroidPoints2dDouble");
    centroid_points_wrapper<double, vec3<double> >("CentroidPoints3dDouble");

    centroid_image_2d_wrapper<float>("CentroidImage2dFloat");
    centroid_image_3d_wrapper<float>("CentroidImage3dFloat");
    centroid_image_2d_wrapper<double>("CentroidImage2dDouble");
    centroid_image_3d_wrapper<double>("CentroidImage3dDouble");

    centroid_masked_image_2d_wrapper<float>("CentroidMaskedImage2dFloat");
    centroid_masked_image_3d_wrapper<float>("CentroidMaskedImage3dFloat");
    centroid_masked_image_2d_wrapper<double>("CentroidMaskedImage2dDouble");
    centroid_masked_image_3d_wrapper<double>("CentroidMaskedImage3dDouble");

    centroid_suite<float>();
    centroid_suite<double>();

    def("centroid_bias_sq", &centroid_bias_sq);
  }

}}}  // namespace dials::algorithms::boost_python
