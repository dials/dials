#include <boost/python.hpp>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/flex_types.h>
#include <cctype>
#include <dxtbx/model/panel.h>
#include <iostream>

namespace recviewer { namespace ext {
  using namespace scitbx;

  typedef scitbx::af::flex<vec3<double> >::type flex_vec3_double;
  typedef scitbx::af::flex<vec2<double> >::type flex_vec2_double;

  static af::shared<vec2<double> > get_target_pixels(dxtbx::model::Panel panel,
                                                     vec3<double> s0,
                                                     int xlim,
                                                     int ylim,
                                                     double maxres) {
    af::shared<vec2<double> > ret;
    vec2<double> xy;

    for (size_t x = 0; x < xlim; x++) {
      for (size_t y = 0; y < ylim; y++) {
        xy[0] = x;
        xy[1] = y;
        try {
          if (panel.get_resolution_at_pixel(s0, xy) > maxres) {
            ret.push_back(xy);
          }
        } catch (...) {
          // FIXME: How can I avoid this? Want to ignore beam center where res = INF.
        }
      }
    }
    return ret;
  }

  static void fill_voxels(const af::flex_int &image,
                          af::flex_double &grid,
                          af::flex_int &cnts,
                          const flex_vec3_double &rotated_S,
                          const flex_vec2_double &xy,
                          const double rec_range) {
    int npoints = grid.accessor().all()[0];  // is this OK?
    double step = 2 * rec_range / npoints;

    for (int i = 0, ilim = xy.size(); i < ilim; i++) {
      int ind_x = rotated_S[i][0] / step + npoints / 2 + 0.5;
      int ind_y = rotated_S[i][1] / step + npoints / 2 + 0.5;
      int ind_z = rotated_S[i][2] / step + npoints / 2 + 0.5;
      int x = xy[i][0];
      int y = xy[i][1];

      if (ind_x >= npoints || ind_y >= npoints || ind_z >= npoints || ind_x < 0
          || ind_y < 0 || ind_z < 0)
        continue;
      grid(ind_x, ind_y, ind_z) += image(y, x);
      cnts(ind_x, ind_y, ind_z)++;
    }
  }

  static void normalize_voxels(af::flex_double &grid, af::flex_int &cnts) {
    for (int i = 0, ilim = grid.size(); i < ilim; i++) {
      if (cnts[i] != 0) {
        grid[i] /= cnts[i];
      }
    }
  }

  void init_module() {
    using namespace boost::python;
    def("get_target_pixels", get_target_pixels);
    def("fill_voxels", fill_voxels);
    def("normalize_voxels", normalize_voxels);
  }

}}  // namespace recviewer::ext

BOOST_PYTHON_MODULE(recviewer_ext) {
  recviewer::ext::init_module();
}
