#ifndef DIALS_ALGORITHMS_IMAGE_DISTORTION_ELLIPSE_H
#define DIALS_ALGORITHMS_IMAGE_DISTORTION_ELLIPSE_H

#include <scitbx/array_family/ref.h>
#include <scitbx/array_family/shared.h>

#include <scitbx/array_family/accessors/c_grid.h>
#include <scitbx/array_family/versa.h>
#include <dxtbx/model/detector.h>

namespace dials { namespace algorithms {

  using dxtbx::model::Panel;

  class CreateEllipticalDistortionMaps {
  public:
    CreateEllipticalDistortionMaps(const Panel &panel) : panel_(panel) {
      std::size_t xsize = panel_.get_image_size()[0];
      std::size_t ysize = panel_.get_image_size()[1];
      dx_.resize(scitbx::af::c_grid<2>(ysize, xsize));
      dy_.resize(scitbx::af::c_grid<2>(ysize, xsize));
    };

    scitbx::af::versa<double, scitbx::af::c_grid<2>> get_dx() const {
      return dx_;
    }

    scitbx::af::versa<double, scitbx::af::c_grid<2>> get_dy() const {
      return dy_;
    }

  private:
    const Panel &panel_;

    scitbx::af::versa<double, scitbx::af::c_grid<2>> dx_;
    scitbx::af::versa<double, scitbx::af::c_grid<2>> dy_;
  };

}}  // namespace dials::algorithms

#endif  // DIALS_ALGORITHMS_IMAGE_DISTORTION_ELLIPSE_H