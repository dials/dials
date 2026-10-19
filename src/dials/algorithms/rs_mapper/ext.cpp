#include <boost/python.hpp>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/flex_types.h>
#include <cctype>
#include <dxtbx/model/panel.h>
#include <iostream>
#include <dials/util/thread_pool.h>
#include <dxtbx/model/experiment_list.h>
#include <dxtbx/model/experiment.h>
#include <dxtbx/model/beam.h>
#include <dxtbx/model/detector.h>
#include <dxtbx/model/scan.h>
#include <dxtbx/imageset.h>
#include <dxtbx/format/image.h>
#include <scitbx/constants.h>

#include <scitbx/vec3.h>
#include <scitbx/vec2.h>

#include <scitbx/mat3.h>
#include <cmath>
#include <dials/error.h>

namespace recviewer { namespace ext {
  using namespace scitbx;
  using dxtbx::format::Image;
  using namespace dxtbx::model;
  using namespace scitbx::af;
  using dxtbx::ImageSequence;
  using scitbx::constants::m_n;
  using scitbx::constants::Planck;
  using namespace boost::python;

  typedef scitbx::af::flex<vec3<double>>::type flex_vec3_double;
  typedef scitbx::af::flex<vec2<double>>::type flex_vec2_double;

  static af::shared<vec2<double>> get_target_pixels(dxtbx::model::Panel panel,
                                                    vec3<double> s0,
                                                    int nfast,
                                                    int nslow,
                                                    double maxres) {
    af::shared<vec2<double>> ret;
    vec2<double> xy;

    for (size_t y = 0; y < nslow; y++) {
      for (size_t x = 0; x < nfast; x++) {
        xy[0] = x;
        xy[1] = y;

        // get_resolution_at_pixel() no longer returns INF, so this is safe
        // Expects coord given in terms of (fast, slow), which is column, row...
        if (panel.get_resolution_at_pixel(s0, xy) > maxres) {
          ret.push_back(xy);
        }
      }
    }
    return ret;
  }

  template <typename ImageType>
  static void fill_voxels(const ImageType& image,
                          af::flex_double& grid,
                          af::flex_int& counts,
                          const flex_vec3_double& rotated_S,
                          const flex_vec2_double& xy,
                          const double rec_range) {
    int npoints = grid.accessor().all()[0];
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
      counts(ind_x, ind_y, ind_z)++;
    }
  }

  static void normalize_voxels(af::flex_double& grid, af::flex_int& counts) {
    for (int i = 0, ilim = grid.size(); i < ilim; i++) {
      if (counts[i] != 0) {
        grid[i] /= counts[i];
      }
    }
  }

  void process_tof_experiment_list(const ExperimentList& experiments,
                                   const double max_res,
                                   af::flex_double& grid,
                                   af::flex_int& counts,
                                   int n_threads,
                                   bool record_all_counts = false) {
    /**
     * Updates grid and counts for all pixels in experiments
     * This is distinct from fill_voxels as the resolution and
     * reciprocal position must be done for each pixel for ToF data
     *
     * If record_all_counts, counts are accepted regardless of max_res
     */

    // Grid geometry
    int npoints = grid.accessor().all()[0];
    double rec_range = 1 / (max_res + 1e-3);
    double step = 2.0 * rec_range / npoints;
    const int total = npoints * npoints * npoints;

    for (std::size_t iexp = 0; iexp < experiments.size(); ++iexp) {
      const Experiment& experiment = experiments[iexp];
      DIALS_ASSERT(experiment.get_type() == TOF);

      // Scan params
      const Scan scan = *experiment.get_scan();
      af::shared<double> img_tof = scan.get_property<double>("time_of_flight");

      // Beam params
      std::shared_ptr<dxtbx::model::BeamBase> beam_ptr = experiment.get_beam();
      std::shared_ptr<PolychromaticBeam> beam =
        std::dynamic_pointer_cast<PolychromaticBeam>(beam_ptr);
      DIALS_ASSERT(beam != nullptr);

      vec3<double> unit_s0 = beam->get_unit_s0();
      double sample_to_source_distance = beam->get_sample_to_source_distance();

      // Detector, Imageset params
      object imageset_obj = *experiment.get_imageset();
      ImageSequence imageset = extract<ImageSequence>(imageset_obj);
      const Detector& detector = *experiment.get_detector();
      std::size_t nframes = imageset.size();

      bool has_goniometer = (experiment.get_goniometer() != nullptr);
      mat3<double> inv_setting_rotation;
      if (has_goniometer) {
        inv_setting_rotation =
          mat3<double>(experiment.get_goniometer()->get_setting_rotation()).inverse();
      }

      // Per-thread grid/counts: shape (n_threads, total), row-major
      versa<double, c_grid<2>> tl_grid(c_grid<2>(n_threads, total), 0.0);
      versa<int, c_grid<2>> tl_counts(c_grid<2>(n_threads, total), 0);

      dials::util::ThreadPool pool(n_threads);

      for (std::size_t frame = 0; frame < nframes; ++frame) {
        double tof = img_tof[frame] * 1e-6;  // (s)
        Image<double> image_data = imageset.get_corrected_data(frame);

        DIALS_ASSERT(image_data.n_tiles() == detector.size());

        for (std::size_t ipanel = 0; ipanel < detector.size(); ++ipanel) {
          const Panel& panel = detector[ipanel];
          auto image = image_data.tile(ipanel).data();
          int nslow = image.accessor()[0];
          int nfast = image.accessor()[1];

          // Distribute y rows across threads
          std::size_t chunk_size = ((std::size_t)nslow + n_threads - 1) / n_threads;
          for (int t = 0; t < n_threads; ++t) {
            std::size_t y_start = t * chunk_size;
            std::size_t y_end = std::min(y_start + chunk_size, (std::size_t)nslow);
            if (y_start >= (std::size_t)nslow) break;

            pool.post([&, tof, nfast, y_start, y_end, t]() {
              for (int y = y_start; y < (int)y_end; ++y) {
                for (int x = 0; x < nfast; ++x) {
                  vec2<double> xy(x, y);
                  vec3<double> s1 = panel.get_pixel_lab_coord(xy);

                  double distance =
                    (s1.length() + sample_to_source_distance) * 1e-3;    // (m)
                  double wl = (Planck * tof) / (m_n * distance) * 1e10;  // (Å)

                  vec3<double> s0 = unit_s0 * (1.0 / wl);
                  s1 = s1 / s1.length() * (1.0 / wl);
                  vec3<double> S = s1 - s0;

                  if (has_goniometer) {
                    S = inv_setting_rotation * S;
                  }

                  int ind_x = static_cast<int>(S[0] / step + npoints / 2 + 0.5);
                  int ind_y = static_cast<int>(S[1] / step + npoints / 2 + 0.5);
                  int ind_z = static_cast<int>(S[2] / step + npoints / 2 + 0.5);

                  if (ind_x < 0 || ind_y < 0 || ind_z < 0 || ind_x >= npoints
                      || ind_y >= npoints || ind_z >= npoints)
                    continue;

                  int idx = ind_x * npoints * npoints + ind_y * npoints + ind_z;
                  if (record_all_counts) {
                    tl_counts(t, idx)++;
                    if (panel.get_resolution_at_pixel(s0, xy) >= max_res) {
                      tl_grid(t, idx) += image(y, x);
                    }
                  } else if (panel.get_resolution_at_pixel(s0, xy) >= max_res) {
                    tl_grid(t, idx) += image(y, x);
                    tl_counts(t, idx)++;
                  }
                }
              }
            });
          }
          pool.wait();
        }
      }

      // Merge threads
      for (int i = 0; i < total; ++i) {
        double g = 0.0;
        int c = 0;
        for (int t = 0; t < n_threads; ++t) {
          g += tl_grid(t, i);
          c += tl_counts(t, i);
        }
        grid[i] += g;
        counts[i] += c;
      }
    }
  }

  // Define specific pointers to the template instances
  auto fill_voxels_int = &fill_voxels<af::flex_int>;
  auto fill_voxels_double = &fill_voxels<af::flex_double>;

  void init_module() {
    using namespace boost::python;
    def("get_target_pixels", get_target_pixels);
    def("fill_voxels", fill_voxels_int);
    def("fill_voxels", fill_voxels_double);
    def("normalize_voxels", normalize_voxels);
    def("process_tof_experiment_list", process_tof_experiment_list);
  }

}}  // namespace recviewer::ext

BOOST_PYTHON_MODULE(recviewer_ext) {
  recviewer::ext::init_module();
}
