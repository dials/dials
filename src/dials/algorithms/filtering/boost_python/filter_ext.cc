/*
 * filter_ext.cc
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
#include <dials/algorithms/filtering/filter.h>

namespace dials { namespace algorithms { namespace filter { namespace boost_python {

  using namespace boost::python;

  using scitbx::af::flex_bool;

  af::shared<bool> by_detector_mask_multipanel_wrapper(
    const af::const_ref<std::size_t> &panel,
    const af::const_ref<int6> bboxes,
    boost::python::tuple mask_tuple,
    int2 scan_range) {
    af::shared<af::const_ref<bool, af::c_grid<2> > > mask(len(mask_tuple));
    for (std::size_t i = 0; i < mask.size(); ++i) {
      flex_bool temp = extract<flex_bool>(mask_tuple[i]);
      DIALS_ASSERT(temp.accessor().all().size() == 2);
      mask[i] = af::const_ref<bool, af::c_grid<2> >(
        &temp[0], af::c_grid<2>(temp.accessor().all()[0], temp.accessor().all()[1]));
    }

    return by_detector_mask_multipanel(panel, bboxes, mask.const_ref(), scan_range);
  }

  void export_is_zeta_valid() {
    def("is_zeta_valid",
        (bool (*)(vec3<double>, vec3<double>, vec3<double>, double)) & is_zeta_valid,
        (arg("m2"), arg("s0"), arg("s1"), arg("zeta_min")));
    def("is_zeta_valid",
        (bool (*)(const CoordinateSystem &, double)) & is_zeta_valid,
        (arg("cs"), arg("zeta_min")));
    def("is_zeta_valid",
        (bool (*)(const Goniometer &, const BeamBase &, vec3<double>, double))
          & is_zeta_valid,
        (arg("g"), arg("b"), arg("s1"), arg("zeta_min")));
  }

  void export_is_xds_small_angle_valid() {
    def("is_xds_small_angle_valid",
        (bool (*)(vec3<double>, vec3<double>, vec3<double>, double))
          & is_xds_small_angle_valid,
        (arg("m2"), arg("s0"), arg("s1"), arg("delta_m")));
    def("is_xds_small_angle_valid",
        (bool (*)(const CoordinateSystem &, double)) & is_xds_small_angle_valid,
        (arg("cs"), arg("delta_m")));
    def("is_xds_small_angle_valid",
        (bool (*)(const Goniometer &, const BeamBase &, vec3<double>, double))
          & is_xds_small_angle_valid,
        (arg("g"), arg("b"), arg("s1"), arg("delta_m")));
  }

  void export_is_xds_angle_valid() {
    def(
      "is_xds_angle_valid",
      (bool (*)(vec3<double>, vec3<double>, vec3<double>, double)) & is_xds_angle_valid,
      (arg("m2"), arg("s0"), arg("s1"), arg("delta_m")));
    def("is_xds_angle_valid",
        (bool (*)(const CoordinateSystem &, double)) & is_xds_angle_valid,
        (arg("cs"), arg("delta_m")));
    def("is_xds_angle_valid",
        (bool (*)(const Goniometer &, const BeamBase &, vec3<double>, double))
          & is_xds_angle_valid,
        (arg("g"), arg("b"), arg("s1"), arg("delta_m")));
  }

  void export_filter_list() {
    def("by_zeta", &by_zeta, (arg("g"), arg("b"), arg("r"), arg("min_zeta")));
    def("by_xds_small_angle",
        &by_xds_small_angle,
        (arg("g"), arg("b"), arg("r"), arg("delta_m")));
    def("by_xds_angle", &by_xds_angle, (arg("g"), arg("b"), arg("r"), arg("delta_m")));
    def(
      "by_bbox_volume",
      (af::shared<bool>(*)(const af::const_ref<int6> &, std::size_t)) & by_bbox_volume,
      (arg("bbox"), arg("num")));
    def("by_bbox_volume",
        (af::shared<bool>(*)(const af::const_ref<int6> &)) & by_bbox_volume,
        (arg("bbox")));

    def("is_bbox_outside_image_range", &is_bbox_outside_image_range);
    def("does_bbox_contain_bad_pixels", &does_bbox_contain_bad_pixels);
    def("is_bbox_valid", &is_bbox_valid);
    def("by_detector_mask", &by_detector_mask);
    def("by_detector_mask", &by_detector_mask_multipanel_wrapper);

    def("by_resolution_at_centroid", &by_resolution_at_centroid);

    def("by_shoebox_mask", &by_shoebox_mask);
  }

  BOOST_PYTHON_MODULE(dials_algorithms_filter_ext) {
    export_is_zeta_valid();
    export_is_xds_small_angle_valid();
    export_is_xds_angle_valid();
    export_filter_list();
  }

}}}}  // namespace dials::algorithms::filter::boost_python
