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
#include <dials/model/data/image_volume.h>
#include <dials/array_family/reflection_table.h>

namespace dials { namespace algorithms {

  using dials::model::Centroid;
  using dials::model::ImageVolume;
  using dials::model::MultiPanelImageVolume;
  using dials::model::Shoebox;
  using dxtbx::model::Detector;
  using dxtbx::model::Panel;
  using dxtbx::model::Scan;

  /**
   * Compute the centroid from a single reflection
   */
  template <typename FloatType>
  Centroid centroid_image_volume(std::size_t index,
                                 int6 bbox,
                                 ImageVolume<FloatType> volume) {
    typedef CentroidMaskedImage3d<FloatType> Centroider;

    // The mask code to use
    int mask_code = Valid | Foreground;

    // Trim the bbox
    bbox = volume.trim_bbox(bbox);

    // Get some arrays
    af::versa<FloatType, af::c_grid<3> > data = volume.extract_data(bbox);
    af::versa<FloatType, af::c_grid<3> > bgrd = volume.extract_background(bbox);
    af::versa<int, af::c_grid<3> > mask = volume.extract_mask(bbox, index);

    // Compute the foreground boolean mask and background substracted data
    af::versa<FloatType, af::c_grid<3> > foreground_data(mask.accessor());
    af::versa<bool, af::c_grid<3> > foreground_mask(mask.accessor());
    for (std::size_t i = 0; i < mask.size(); ++i) {
      foreground_data[i] = data[i] - bgrd[i];
      foreground_mask[i] =
        ((mask[i] & mask_code) == mask_code) && (foreground_data[i] >= 0);
    }

    // Set the result
    Centroid result;
    try {
      Centroider algorithm(foreground_data.const_ref(), foreground_mask.const_ref());
      result.px.position = algorithm.mean() + vec3<double>(bbox[0], bbox[2], bbox[4]);
      result.px.variance = algorithm.variance();
      result.px.std_err_sq = algorithm.mean_sq_error();
    } catch (dials::error) {
      double xmid = (bbox[1] + bbox[0]) / 2.0;
      double ymid = (bbox[3] + bbox[2]) / 2.0;
      double zmid = (bbox[5] + bbox[4]) / 2.0;
      result.px.position = vec3<double>(xmid, ymid, zmid);
      result.px.variance = vec3<double>(0, 0, 0);
      result.px.std_err_sq = vec3<double>(1.0 / 12.0, 1.0 / 12.0, 1.0 / 12.0);
    }
    return result;
  }

  class Centroider {
  public:
    Centroider() {}

    void add(const Detector& detector) {
      detector_.push_back(detector);
      scan_.push_back(boost::none);
    }

    void add(const Detector& detector, const Scan& scan) {
      detector_.push_back(detector);
      scan_.push_back(scan);
    }

    void shoebox(af::reflection_table reflections) const {
      // Check stuff
      DIALS_ASSERT(reflections.is_consistent());
      DIALS_ASSERT(reflections.size() > 0);
      DIALS_ASSERT(reflections.contains("shoebox"));
      DIALS_ASSERT(reflections.contains("id"));
      DIALS_ASSERT(detector_.size() > 0);
      DIALS_ASSERT(detector_.size() == scan_.size());

      // Loop through all the reflections
      af::const_ref<int> id = reflections["id"];
      af::const_ref<Shoebox<> > shoebox = reflections["shoebox"];
      af::ref<vec3<double> > xyzobs_px_value = reflections["xyzobs.px.value"];
      af::ref<vec3<double> > xyzobs_px_variance = reflections["xyzobs.px.variance"];
      af::ref<vec3<double> > xyzobs_mm_value = reflections["xyzobs.mm.value"];
      af::ref<vec3<double> > xyzobs_mm_variance = reflections["xyzobs.mm.variance"];
      for (std::size_t i = 0; i < shoebox.size(); ++i) {
        // Compute the pixel centroid
        Centroid centroid = shoebox[i].centroid_foreground_minus_background();

        // Get the panel
        DIALS_ASSERT(id[i] >= 0);
        DIALS_ASSERT(id[i] < detector_.size());
        const Detector& d = detector_[id[i]];
        DIALS_ASSERT(shoebox[i].panel < d.size());
        const Panel& p = d[shoebox[i].panel];

        // Get the mm centroid
        vec2<double> mm = p.pixel_to_millimeter(
          vec2<double>(centroid.px.position[0], centroid.px.position[1]));
        vec2<double> pixel_size = p.get_pixel_size();
        double xscale = pixel_size[0] * pixel_size[0];
        double yscale = pixel_size[1] * pixel_size[1];
        centroid.mm.position[0] = mm[0];
        centroid.mm.position[1] = mm[1];
        centroid.mm.std_err_sq[0] = centroid.px.std_err_sq[0] * xscale;
        centroid.mm.std_err_sq[1] = centroid.px.std_err_sq[1] * yscale;

        // Get the phi centroid
        double zscale = 0.0;
        double phi = 0.0;
        if (scan_[id[i]]) {
          phi = scan_[id[i]]->get_angle_from_array_index(centroid.px.position[2]);
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

    template <typename FloatType>
    void volume(af::reflection_table reflections,
                MultiPanelImageVolume<FloatType> volume) const {
      // Check stuff
      DIALS_ASSERT(reflections.is_consistent());
      DIALS_ASSERT(reflections.size() > 0);
      DIALS_ASSERT(reflections.contains("bbox"));
      DIALS_ASSERT(reflections.contains("id"));
      DIALS_ASSERT(reflections.contains("panel"));
      DIALS_ASSERT(detector_.size() > 0);
      DIALS_ASSERT(detector_.size() == scan_.size());

      // Loop through all the reflections
      af::const_ref<int> id = reflections["id"];
      af::const_ref<std::size_t> panel = reflections["panel"];
      af::const_ref<int6> bbox = reflections["bbox"];
      af::ref<vec3<double> > xyzobs_px_value = reflections["xyzobs.px.value"];
      af::ref<vec3<double> > xyzobs_px_variance = reflections["xyzobs.px.variance"];
      af::ref<vec3<double> > xyzobs_mm_value = reflections["xyzobs.mm.value"];
      af::ref<vec3<double> > xyzobs_mm_variance = reflections["xyzobs.mm.variance"];
      for (std::size_t i = 0; i < bbox.size(); ++i) {
        // Compute the centroid
        Centroid centroid = centroid_image_volume(i, bbox[i], volume.get(panel[i]));

        // Get the panel
        DIALS_ASSERT(id[i] >= 0);
        DIALS_ASSERT(id[i] < detector_.size());
        const Detector& d = detector_[id[i]];
        DIALS_ASSERT(panel[i] < d.size());
        const Panel& p = d[panel[i]];

        // Get the mm centroid
        vec2<double> mm = p.pixel_to_millimeter(
          vec2<double>(centroid.px.position[0], centroid.px.position[1]));
        vec2<double> pixel_size = p.get_pixel_size();
        double xscale = pixel_size[0] * pixel_size[0];
        double yscale = pixel_size[1] * pixel_size[1];
        centroid.mm.position[0] = mm[0];
        centroid.mm.position[1] = mm[1];
        centroid.mm.std_err_sq[0] = centroid.px.std_err_sq[0] * xscale;
        centroid.mm.std_err_sq[1] = centroid.px.std_err_sq[1] * yscale;

        // Get the phi centroid
        double zscale = 0.0;
        double phi = 0.0;
        if (scan_[id[i]]) {
          phi = scan_[id[i]]->get_angle_from_array_index(centroid.px.position[2]);
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

}}  // namespace dials::algorithms

#endif  // DIALS_ALGORITHMS_CENTROID_SIMPLE_ALGORITHM_H
