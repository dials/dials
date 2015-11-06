/*
 * algorithm.cc
 *
 *  copyright (c) 2013 diamond light source
 *
 *  author: james parkhurst
 *
 *  this code is distributed under the bsd license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_CENTROID_SIMPLE_ALGORITHM_H
#define DIALS_ALGORITHMS_CENTROID_SIMPLE_ALGORITHM_H

#include <boost/optional.hpp>
#include <dxtbx/model/detector.h>
#include <dxtbx/model/scan.h>
#include <dials/array_family/reflection_table.h>

namespace dials { namespace algorithms {

  using dxtbx::model::Detector;
  using dxtbx::model::Panel;
  using dxtbx::model::Scan;
  using dials::model::Shoebox;
  using dials::model::Centroid;

  class Centroider {
  public:

    Centroider() {}

    void add(const Detector &detector) {
      detector_.push_back(detector);
      scan_.push_back(boost::none);
    }

    void add(const Detector &detector, const Scan& scan) {
      detector_.push_back(detector);
      scan_.push_back(scan);
    }

    void operator()(af::reflection_table reflections) const {
      // Check stuff
      DIALS_ASSERT(reflections.is_consistent());
      DIALS_ASSERT(reflections.size() > 0);
      DIALS_ASSERT(reflections.contains("shoebox"));
      DIALS_ASSERT(reflections.contains("id"));
      DIALS_ASSERT(detector_.size() > 0);
      DIALS_ASSERT(detector_.size() == scan_.size());

      // Loop through all the reflections
      af::const_ref<std::size_t> id = reflections["id"];
      af::const_ref< Shoebox<> > shoebox = reflections["shoebox"];
      af::ref< vec3<double> > xyzobs_px_value = reflections["xyzobs.px.value"];
      af::ref< vec3<double> > xyzobs_px_variance = reflections["xyzobs.px.variance"];
      af::ref< vec3<double> > xyzobs_mm_value = reflections["xyzobs.mm.value"];
      af::ref< vec3<double> > xyzobs_mm_variance = reflections["xyzobs.mm.variance"];
      for (std::size_t i = 0; i < shoebox.size(); ++i) {

        // Compute the pixel centroid
        Centroid centroid = shoebox[i].centroid_foreground_minus_background();

        // Get the panel
        DIALS_ASSERT(id[i] < detector_.size());
        const Detector& d = detector_[id[i]];
        DIALS_ASSERT(shoebox[i].panel < d.size());
        const Panel& p = d[shoebox[i].panel];

        // Get the mm centroid
        vec2<double> mm = p.pixel_to_millimeter(vec2<double>(
              centroid.px.position[0],
              centroid.px.position[1]));
        vec2<double> pixel_size = p.get_pixel_size();
        double xscale = pixel_size[0]*pixel_size[0];
        double yscale = pixel_size[1]*pixel_size[1];
        centroid.mm.position[0] = mm[0];
        centroid.mm.position[1] = mm[1];
        centroid.mm.std_err_sq[0] = centroid.px.std_err_sq[0] * xscale;
        centroid.mm.std_err_sq[1] = centroid.px.std_err_sq[1] * yscale;

        // Get the phi centroid
        double zscale = 0.0;
        double phi = 0.0;
        if (scan_[id[i]]) {
          phi = scan_[id[i]]->get_angle_from_array_index(
              centroid.px.position[2]);
          zscale = scan_[id[i]]->get_oscillation()[1];
          zscale *= zscale;
        }
        centroid.mm.position[2] = phi;
        centroid.mm.std_err_sq[2] = centroid.px.std_err_sq[2] * zscale;

        // Set array values
        xyzobs_px_value[i] = centroid.px.position;
        xyzobs_mm_value[i] = centroid.mm.position;
        xyzobs_px_variance[i] = centroid.px.std_err_sq;
        xyzobs_mm_variance[i] = centroid.mm.std_err_sq;
      }
    }

  private:

    std::vector<Detector> detector_;
    std::vector<boost::optional<Scan> > scan_;
  };

}} // namespace dials::algorithms

#endif // DIALS_ALGORITHMS_CENTROID_SIMPLE_ALGORITHM_H
