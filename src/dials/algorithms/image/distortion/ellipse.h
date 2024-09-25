#ifndef DIALS_ALGORITHMS_IMAGE_DISTORTION_ELLIPSE_H
#define DIALS_ALGORITHMS_IMAGE_DISTORTION_ELLIPSE_H

#include <scitbx/array_family/accessors/c_grid.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/mat2.h>
#include <scitbx/vec2.h>
#include <dxtbx/model/detector.h>

namespace dials { namespace algorithms {

  using dxtbx::model::Panel;
  using scitbx::mat2;
  using scitbx::vec2;
  using scitbx::vec3;

  class CreateEllipticalDistortionMaps {
  public:
    CreateEllipticalDistortionMaps(const Panel &panel,
                                   const mat2<double> ellipse_matrix,
                                   const vec3<double> fast,
                                   const vec3<double> slow,
                                   const vec3<double> mid) {
      std::size_t xsize = panel.get_image_size()[0];
      std::size_t ysize = panel.get_image_size()[1];
      dx_.resize(scitbx::af::c_grid<2>(ysize, xsize));
      dy_.resize(scitbx::af::c_grid<2>(ysize, xsize));

      for (std::size_t j = 0; j < ysize; ++j) {
        for (std::size_t i = 0; i < xsize; ++i) {
          vec3<double> lab = panel.get_pixel_lab_coord(vec2<double>(i + 0.5, j + 0.5));
          vec3<double> offset = lab - mid;
          double x = offset * fast;  // undistorted X coordinate (mm)
          double y = offset * slow;  // undistorted Y coordinate (mm)
          vec2<double> distort =
            ellipse_matrix * vec2<double>(x, y);  // distorted by ellipse matrix

          // store correction in units of the pixel size
          dx_(j, i) = (x - distort[0]) / panel.get_pixel_size()[0];
          dy_(j, i) = (y - distort[1]) / panel.get_pixel_size()[1];
        }
      }
    };

    scitbx::af::versa<double, scitbx::af::c_grid<2>> get_dx() const {
      return dx_;
    }

    scitbx::af::versa<double, scitbx::af::c_grid<2>> get_dy() const {
      return dy_;
    }

  private:
    scitbx::af::versa<double, scitbx::af::c_grid<2>> dx_;
    scitbx::af::versa<double, scitbx::af::c_grid<2>> dy_;
  };

}}  // namespace dials::algorithms

#endif  // DIALS_ALGORITHMS_IMAGE_DISTORTION_ELLIPSE_H