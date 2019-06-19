#ifndef DIALS_UTIL_MASKING_GONIOMETER_SHADOW_MASKING_H
#define DIALS_UTIL_MASKING_GONIOMETER_SHADOW_MASKING_H

#include <algorithm>
#include <boost/geometry.hpp>
#include <boost/geometry/algorithms/assign.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/adapted/boost_tuple.hpp>
#include <scitbx/math/r3_rotation.h>
#include <dxtbx/format/image.h>
#include <dxtbx/model/detector.h>
#include <dxtbx/model/panel.h>
#include <dxtbx/model/multi_axis_goniometer.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/util/masking.h>
#include <dials/error.h>

BOOST_GEOMETRY_REGISTER_BOOST_TUPLE_CS(boost::geometry::cs::cartesian)

namespace dials { namespace util { namespace masking {

  using dxtbx::format::Image;
  using dxtbx::format::ImageTile;
  using dxtbx::model::Detector;
  using dxtbx::model::MultiAxisGoniometer;
  using dxtbx::model::Panel;
  using scitbx::vec2;
  using scitbx::vec3;

  /**
   * A class to mask multiple resolution ranges
   */
  class GoniometerShadowMaskGenerator {
  public:
    /**
     * Initialise the resolution at each pixel
     * @param beam The beam model
     * @param panel The panel model
     */
    GoniometerShadowMaskGenerator(
      const MultiAxisGoniometer &goniometer,
      const scitbx::af::const_ref<vec3<double> > &extrema_at_datum,
      const scitbx::af::const_ref<std::size_t> &axis)
        : goniometer_(goniometer),
          extrema_at_datum_(extrema_at_datum.begin(), extrema_at_datum.end()),
          axis_(axis.begin(), axis.end()) {}

    scitbx::af::shared<vec3<double> > extrema_at_scan_angle(double scan_angle) const {
      scitbx::af::shared<vec3<double> > axes = goniometer_.get_axes();
      scitbx::af::shared<double> angles = goniometer_.get_angles();
      std::size_t scan_axis = goniometer_.get_scan_axis();
      angles[scan_axis] = scan_angle;
      scitbx::af::shared<vec3<double> > extrema(extrema_at_datum_.begin(),
                                                extrema_at_datum_.end());

      for (std::size_t i = 0; i < axes.size(); i++) {
        scitbx::mat3<double> rotation =
          scitbx::math::r3_rotation::axis_and_angle_as_matrix(axes[i], angles[i], true);
        for (std::size_t j = 0; j < axis_.size(); j++) {
          if (axis_[j] > i) {
            continue;
          }
          extrema[j] = rotation * extrema[j];
        }
      }
      return extrema;
    }

    scitbx::af::shared<scitbx::af::shared<vec2<double> > > project_extrema(
      const Detector &detector,
      double scan_angle) const {
      scitbx::af::shared<vec3<double> > coords = extrema_at_scan_angle(scan_angle);

      typedef boost::tuple<double, double> point_t;
      typedef boost::geometry::model::polygon<point_t> polygon_t;

      scitbx::af::shared<scitbx::af::shared<vec2<double> > > result;

      std::vector<point_t> points;
      for (std::size_t i = 0; i < detector.size(); i++) {
        Panel panel = detector[i];

        /* project coordinates onto panel plane */
        for (std::size_t j = 0; j < coords.size(); j++) {
          vec3<double> coord = panel.get_D_matrix() * coords[j];
          double z = coord[2];
          if (z > 0) {
            point_t p(coord[0] / z, coord[1] / z);
            points.push_back(p);
          }
        }
        polygon_t poly;
        boost::geometry::assign_points(poly, points);

        scitbx::af::shared<vec2<double> > shadow_points;

        if (points.size() == 0) {
          result.push_back(shadow_points);
          continue;
        }

        // Construct detector polygon - points should be clockwise
        std::vector<point_t> corners;
        corners.push_back(point_t(0, 0));
        corners.push_back(point_t(0, panel.get_image_size_mm()[1]));
        corners.push_back(
          point_t(panel.get_image_size_mm()[0], panel.get_image_size_mm()[1]));
        corners.push_back(point_t(panel.get_image_size_mm()[0], 0));
        corners.push_back(point_t(0, 0));
        polygon_t det;
        boost::geometry::assign_points(det, corners);

        // Check the validity of the polygon
        boost::geometry::validity_failure_type failure;
        bool valid = boost::geometry::is_valid(poly, failure);

        // if the invalidity is only due to lack of closing points and/or wrongly
        // oriented rings, then bg::correct can fix it
        bool could_be_fixed = (failure == boost::geometry::failure_not_closed
                               || boost::geometry::failure_wrong_orientation);
        if (!valid) {
          if (could_be_fixed) {
            boost::geometry::correct(poly);
            valid = boost::geometry::is_valid(poly);
          }
        }
        if (!valid) {
          std::cout << "Invalid polygon geometry: " << boost::geometry::dsv(poly)
                    << std::endl;
        }
        DIALS_ASSERT(valid);

        polygon_t shadow;
        boost::geometry::convex_hull(poly, shadow);

        // Compute the intersection of the shadow in the detector plane with the
        // detector
        std::deque<polygon_t> output;
        boost::geometry::intersection(det, shadow, output);

        // Extract the coordinates of the shadow on the detector, and convert from mm
        // to pixel coordinates
        if (output.size()) {
          vec2<double> px = panel.get_pixel_size();
          polygon_t hull = output[0];
          for (auto shadow_iterator = hull.outer().begin();
               shadow_iterator != hull.outer().end();
               ++shadow_iterator) {
            vec2<double> p(boost::geometry::get<0>(*shadow_iterator) / px[0],
                           boost::geometry::get<1>(*shadow_iterator) / px[1]);
            shadow_points.push_back(p);
          }
        }
        result.push_back(shadow_points);
      }
      return result;
    }

    Image<bool> get_mask(const Detector &detector, double scan_angle) const {
      scitbx::af::shared<scitbx::af::shared<vec2<double> > > shadow_boundary =
        project_extrema(detector, scan_angle);

      Image<bool> mask;

      for (std::size_t i = 0; i < detector.size(); i++) {
        Panel panel = detector[i];

        vec2<std::size_t> image_size = panel.get_image_size();
        scitbx::af::versa<bool, scitbx::af::c_grid<2> > mask_data(
          scitbx::af::c_grid<2>(image_size[1], image_size[0]), true);

        ImageTile<bool> mask_tile(mask_data);

        if (shadow_boundary[i].size() >= 3) {
          dials::util::mask_untrusted_polygon(mask_data.ref(),
                                              shadow_boundary[i].const_ref());
        }
        mask.push_back(mask_tile);
      }

      return mask;
    }

  private:
    const MultiAxisGoniometer &goniometer_;
    scitbx::af::shared<vec3<double> > extrema_at_datum_;
    scitbx::af::shared<std::size_t> axis_;
  };

}}}  // namespace dials::util::masking

#endif /* DIALS_UTIL_MASKING_GONIOMETER_SHADOW_MASKING_H */
