/*
 * modeller.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */

#ifndef DIALS_ALGORITHMS_PROFILE_MODEL_GAUSSIAN_RS_MODELLER_H
#define DIALS_ALGORITHMS_PROFILE_MODEL_GAUSSIAN_RS_MODELLER_H

#include <algorithm>
#include <cmath>
#include <memory>
#include <fstream>
#include <dials/algorithms/profile_model/gaussian_rs/transform/transform.h>
#include <dials/algorithms/profile_model/modeller/empirical_modeller.h>
#include <dials/algorithms/profile_model/modeller/single_sampler.h>
#include <dials/algorithms/profile_model/modeller/grid_sampler.h>
#include <dials/algorithms/profile_model/modeller/circle_sampler.h>
#include <dials/algorithms/profile_model/modeller/ewald_sphere_sampler.h>
#include <dials/algorithms/integration/fit/fitting.h>

namespace dials { namespace algorithms {

  using dials::algorithms::profile_model::gaussian_rs::CoordinateSystem;
  using dials::algorithms::profile_model::gaussian_rs::transform::TransformForward;
  using dials::algorithms::profile_model::gaussian_rs::transform::TransformReverse;
  using dials::algorithms::profile_model::gaussian_rs::transform::TransformSpec;
  using dials::model::Shoebox;
  using dxtbx::model::BeamBase;
  using dxtbx::model::Detector;
  using dxtbx::model::Goniometer;
  using dxtbx::model::Scan;

  /**
   * A base class to initialize the sampler
   */
  class GaussianRSProfileModellerBase {
  public:
    enum GridMethod {
      Single = 1,
      RegularGrid = 2,
      CircularGrid = 3,
      SphericalGrid = 4,
    };

    enum FitMethod { ReciprocalSpace = 1, DetectorSpace = 2, CellCache = 3 };

    GaussianRSProfileModellerBase(const std::shared_ptr<BeamBase> beam,
                                  const Detector& detector,
                                  const Goniometer& goniometer,
                                  const Scan& scan,
                                  double sigma_b,
                                  double sigma_m,
                                  double n_sigma,
                                  std::size_t grid_size,
                                  std::size_t num_scan_points,
                                  int grid_method,
                                  int fit_method)
        : beam_(beam),
          detector_(detector),
          goniometer_(goniometer),
          scan_(scan),
          sigma_b_(sigma_b),
          sigma_m_(sigma_m),
          n_sigma_(n_sigma),
          grid_size_(grid_size),
          num_scan_points_(num_scan_points),
          grid_method_(grid_method),
          fit_method_(fit_method),
          sampler_(init_sampler(beam,
                                detector,
                                goniometer,
                                scan,
                                num_scan_points,
                                grid_method)) {}

  protected:
    std::shared_ptr<SamplerIface> init_sampler(std::shared_ptr<BeamBase> beam,
                                               const Detector& detector,
                                               const Goniometer& goniometer,
                                               const Scan& scan,
                                               std::size_t num_scan_points,
                                               int grid_method) {
      int2 scan_range = scan.get_array_range();
      std::shared_ptr<SamplerIface> sampler;
      if (grid_method == RegularGrid || grid_method == CircularGrid) {
        if (detector.size() > 1) {
          grid_method = Single;
        }
      }
      switch (grid_method) {
      case Single:
        sampler = std::make_shared<SingleSampler>(scan_range, num_scan_points);
        break;
      case RegularGrid:
        DIALS_ASSERT(detector.size() == 1);
        sampler = std::make_shared<GridSampler>(
          detector[0].get_image_size(), scan_range, int3(3, 3, num_scan_points));
        break;
      case CircularGrid:
        DIALS_ASSERT(detector.size() == 1);
        sampler = std::make_shared<CircleSampler>(
          detector[0].get_image_size(), scan_range, num_scan_points);
        break;
      case SphericalGrid:
        sampler = std::make_shared<EwaldSphereSampler>(
          beam, detector, goniometer, scan, num_scan_points);
      default:
        throw DIALS_ERROR("Unknown grid method");
      };
      return sampler;
    }

    std::shared_ptr<BeamBase> beam_;
    Detector detector_;
    Goniometer goniometer_;
    Scan scan_;
    double sigma_b_;
    double sigma_m_;
    double n_sigma_;
    std::size_t grid_size_;
    std::size_t num_scan_points_;
    int grid_method_;
    int fit_method_;
    std::shared_ptr<SamplerIface> sampler_;
  };

  namespace detail {

    struct check_mask_code {
      int mask_code;
      check_mask_code(int code) : mask_code(code) {}
      bool operator()(int a) const {
        return ((a & mask_code) == mask_code);
      }
    };

    struct check_either_mask_code {
      int mask_code1;
      int mask_code2;
      check_either_mask_code(int code1, int code2)
          : mask_code1(code1), mask_code2(code2) {}
      bool operator()(int a) const {
        return ((a & mask_code1) == mask_code1) || ((a & mask_code2) == mask_code2);
      }
    };

  }  // namespace detail

  /**
   * The profile modeller for the gaussian rs profile model
   */
  class GaussianRSProfileModeller : public GaussianRSProfileModellerBase,
                                    public EmpiricalProfileModeller {
  public:
    /**
     * Initialize
     * @param beam The beam model
     * @param detector The detector model
     * @param goniometer The goniometer model
     * @param scan The scan model
     * @param sigma_b The beam divergence
     * @param sigma_m The mosaicity
     * @param n_sigma The extent
     * @param grid_size The size of the profile grid
     * @param num_scan_points The number of phi scan points
     * @param threshold The modelling threshold value
     * @param grid_method The gridding method
     */
    GaussianRSProfileModeller(std::shared_ptr<BeamBase> beam,
                              const Detector& detector,
                              const Goniometer& goniometer,
                              const Scan& scan,
                              double sigma_b,
                              double sigma_m,
                              double n_sigma,
                              std::size_t grid_size,
                              std::size_t num_scan_points,
                              double threshold,
                              int grid_method,
                              int fit_method)
        : GaussianRSProfileModellerBase(beam,
                                        detector,
                                        goniometer,
                                        scan,
                                        sigma_b,
                                        sigma_m,
                                        n_sigma,
                                        grid_size,
                                        num_scan_points,
                                        grid_method,
                                        fit_method),
          EmpiricalProfileModeller(
            sampler_->size(),
            int3(2 * grid_size + 1, 2 * grid_size + 1, 2 * grid_size + 1),
            threshold),
          spec_(beam,
                detector,
                goniometer,
                scan,
                sigma_b,
                sigma_m,
                n_sigma,
                grid_size) {
      DIALS_ASSERT(sampler_ != 0);
      cube_cache_.resize(sampler_->size());
      cube_valid_.resize(sampler_->size(), false);
    }

    std::shared_ptr<BeamBase> beam() const {
      return beam_;
    }

    Detector detector() const {
      return detector_;
    }

    Goniometer goniometer() const {
      return goniometer_;
    }

    Scan scan() const {
      return scan_;
    }

    double sigma_b() const {
      return sigma_b_;
    }

    double sigma_m() const {
      return sigma_m_;
    }

    double n_sigma() const {
      return n_sigma_;
    }

    std::size_t grid_size() const {
      return grid_size_;
    }

    std::size_t num_scan_points() const {
      return num_scan_points_;
    }

    double threshold() const {
      return threshold_;
    }

    int grid_method() const {
      return grid_method_;
    }

    int fit_method() const {
      return fit_method_;
    }

    vec3<double> coord(std::size_t index) const {
      return sampler_->coord(index);
    }

    /**
     * Model the profiles from the reflections
     * @param reflections The reflection list
     */
    void model(af::reflection_table reflections) {
      // Check input is OK
      DIALS_ASSERT(reflections.is_consistent());
      DIALS_ASSERT(reflections.contains("shoebox"));
      DIALS_ASSERT(reflections.contains("flags"));
      DIALS_ASSERT(reflections.contains("partiality"));
      DIALS_ASSERT(reflections.contains("s1"));
      DIALS_ASSERT(reflections.contains("xyzcal.px"));
      DIALS_ASSERT(reflections.contains("xyzcal.mm"));

      // Get some data
      af::const_ref<Shoebox<>> sbox = reflections["shoebox"];
      af::const_ref<double> partiality = reflections["partiality"];
      af::const_ref<vec3<double>> s1 = reflections["s1"];
      af::const_ref<vec3<double>> xyzpx = reflections["xyzcal.px"];
      af::const_ref<vec3<double>> xyzmm = reflections["xyzcal.mm"];
      af::ref<std::size_t> flags = reflections["flags"];

      // Loop through all the reflections and add them to the model
      for (std::size_t i = 0; i < reflections.size(); ++i) {
        DIALS_ASSERT(sbox[i].is_consistent());

        // Check if we want to use this reflection
        if (check1(flags[i], partiality[i], sbox[i])) {
          // Create the coordinate system
          vec3<double> m2 = spec_.goniometer().get_rotation_axis();
          vec3<double> s0 = spec_.beam()->get_s0();
          CoordinateSystem cs(m2, s0, s1[i], xyzmm[i][2]);

          // Create the data array
          af::versa<double, af::c_grid<3>> data(sbox[i].data.accessor());
          std::transform(sbox[i].data.begin(),
                         sbox[i].data.end(),
                         sbox[i].background.begin(),
                         data.begin(),
                         std::minus<double>());

          // Create the mask array
          af::versa<bool, af::c_grid<3>> mask(sbox[i].mask.accessor());
          std::transform(sbox[i].mask.begin(),
                         sbox[i].mask.end(),
                         mask.begin(),
                         detail::check_mask_code(Valid | Foreground));

          // Compute the transform
          TransformForward<double> transform(
            spec_, cs, sbox[i].bbox, sbox[i].panel, data.const_ref(), mask.const_ref());

          // Get the indices and weights of the profiles
          af::shared<std::size_t> indices =
            sampler_->nearest_n(sbox[i].panel, xyzpx[i]);
          af::shared<double> weights(indices.size());
          for (std::size_t j = 0; j < indices.size(); ++j) {
            weights[j] = sampler_->weight(indices[j], sbox[i].panel, xyzpx[i]);
          }

          // Add the profile
          add(
            indices.const_ref(), weights.const_ref(), transform.profile().const_ref());

          // Set the flags
          flags[i] |= af::UsedInModelling;
        }
      }
    }

    /**
     * Return a profile fitter
     * @return The profile fitter class
     */
    af::shared<bool> fit(af::reflection_table reflections) const {
      af::shared<bool> success;
      switch (fit_method_) {
      case ReciprocalSpace:
      case CellCache:
        success = fit_reciprocal_space(reflections);
        break;
      case DetectorSpace:
        success = fit_detector_space(reflections);
        break;
      default:
        throw DIALS_ERROR("Unknown fitting method");
      };
      return success;
    }

    /**
     * Return a profile fitter
     * @return The profile fitter class
     */
    void validate(af::reflection_table reflections) const {
      switch (fit_method_) {
      case ReciprocalSpace:
      case CellCache:
        fit_reciprocal_space(reflections);
        break;
      case DetectorSpace:
        fit_detector_space(reflections);
        break;
      default:
        throw DIALS_ERROR("Unknown fitting method");
      };
    }

    /**
     * Return a profile fitter
     * @return The profile fitter class
     */
    af::shared<bool> fit_reciprocal_space(af::reflection_table reflections) const {
      // Check input is OK
      DIALS_ASSERT(reflections.is_consistent());
      DIALS_ASSERT(reflections.contains("shoebox"));
      DIALS_ASSERT(reflections.contains("flags"));
      DIALS_ASSERT(reflections.contains("partiality"));
      DIALS_ASSERT(reflections.contains("s1"));
      DIALS_ASSERT(reflections.contains("xyzcal.px"));
      DIALS_ASSERT(reflections.contains("xyzcal.mm"));

      // Get some data
      af::const_ref<Shoebox<>> sbox = reflections["shoebox"];
      af::const_ref<vec3<double>> s1 = reflections["s1"];
      af::const_ref<vec3<double>> xyzpx = reflections["xyzcal.px"];
      af::const_ref<vec3<double>> xyzmm = reflections["xyzcal.mm"];
      af::ref<std::size_t> flags = reflections["flags"];
      af::ref<double> intensity_val = reflections["intensity.prf.value"];
      af::ref<double> intensity_var = reflections["intensity.prf.variance"];
      af::ref<double> reference_cor = reflections["profile.correlation"];
      // af::ref<double> reference_rmsd = reflections["profile.rmsd"];

      // Loop through all the reflections and process them
      af::shared<bool> success(reflections.size(), false);
      for (std::size_t i = 0; i < reflections.size(); ++i) {
        DIALS_ASSERT(sbox[i].is_consistent());

        // Set values to bad
        intensity_val[i] = 0.0;
        intensity_var[i] = -1.0;
        reference_cor[i] = 0.0;
        // reference_rmsd[i] = 0.0;
        flags[i] &= ~af::IntegratedPrf;
        bool integrate = !(flags[i] & af::DontIntegrate);

        // Check if we want to use this reflection
        if (integrate) {
          try {
            // Get the sampler cell index
            std::size_t index = sampler_->nearest(sbox[i].panel, xyzpx[i]);

            if (fit_method_ == CellCache) {
              // --- CELL-CACHE PATH ---
              // Skip reflections with no learned profile for this cell
              if (!valid(index)) continue;

              // Lazy-init the cube for this cell
              if (!cube_valid_[index]) {
                cube_cache_[index] = build_cell_cube(index, sbox[i].panel);
                cube_valid_[index] = true;
              }
              const CellCube& cc = cube_cache_[index];

              // Compute sub-pixel offset from cell centroid
              vec3<double> cell_coord = sampler_->coord(index);
              double dx = xyzpx[i][0] - cell_coord[0];
              double dy = xyzpx[i][1] - cell_coord[1];
              double dz = xyzpx[i][2] - cell_coord[2];
              // Clamp to cube margin range
              dx = std::max(-1.0, std::min(1.0, dx));
              dy = std::max(-1.0, std::min(1.0, dy));
              dz = std::max(-1.0, std::min(1.0, dz));

              // Extract reference sub-volume via trilinear interpolation
              ExtractedRef ref = extract_from_cube(cc, dx, dy, dz, sbox[i].bbox);

              // Build combined mask (reference mask AND shoebox foreground mask)
              af::versa<bool, af::c_grid<3>> m(ref.mask.accessor());
              for (std::size_t j = 0; j < m.size(); ++j) {
                m[j] =
                  ref.mask[j]
                  && ((sbox[i].mask[j] & (Valid | Foreground)) == (Valid | Foreground));
              }

              // Convert float shoebox data/background to double for ProfileFitter
              af::versa<double, af::c_grid<3>> data_d(sbox[i].data.accessor());
              af::versa<double, af::c_grid<3>> bg_d(sbox[i].background.accessor());
              for (std::size_t j = 0; j < data_d.size(); ++j) {
                data_d[j] = static_cast<double>(sbox[i].data[j]);
                bg_d[j] = static_cast<double>(sbox[i].background[j]);
              }

              // IRLS fit using raw shoebox data and extracted detector-space reference
              ProfileFitter<double> fit(data_d.const_ref(),
                                        bg_d.const_ref(),
                                        m.const_ref(),
                                        ref.profile.const_ref(),
                                        1e-3,
                                        100);

              // Store results
              intensity_val[i] = fit.intensity()[0];
              intensity_var[i] = fit.variance()[0];
              reference_cor[i] = fit.correlation();
              flags[i] |= af::IntegratedPrf;
              success[i] = true;

            } else {
              // --- EXISTING RECIPROCAL-SPACE PATH ---
              data_const_reference p = data(index).const_ref();
              mask_const_reference mask1 = mask(index).const_ref();

              // Create the coordinate system
              vec3<double> m2 = spec_.goniometer().get_rotation_axis();
              vec3<double> s0 = spec_.beam()->get_s0();
              CoordinateSystem cs(m2, s0, s1[i], xyzmm[i][2]);

              // Create the data array
              af::versa<double, af::c_grid<3>> data(sbox[i].data.accessor());
              std::copy(sbox[i].data.begin(), sbox[i].data.end(), data.begin());

              // Create the background array
              af::versa<double, af::c_grid<3>> background(
                sbox[i].background.accessor());
              std::copy(sbox[i].background.begin(),
                        sbox[i].background.end(),
                        background.begin());

              // Create the mask array
              af::versa<bool, af::c_grid<3>> mask(sbox[i].mask.accessor());
              std::transform(sbox[i].mask.begin(),
                             sbox[i].mask.end(),
                             mask.begin(),
                             detail::check_mask_code(Valid | Foreground));

              // Compute the transform
              TransformForward<double> transform(spec_,
                                                 cs,
                                                 sbox[i].bbox,
                                                 sbox[i].panel,
                                                 data.const_ref(),
                                                 background.const_ref(),
                                                 mask.const_ref());

              // Get the transformed shoebox
              data_const_reference c = transform.profile().const_ref();
              data_const_reference b = transform.background().const_ref();
              mask_const_reference mask2 = transform.mask().const_ref();
              af::versa<bool, af::c_grid<3>> m(mask2.accessor());
              DIALS_ASSERT(mask1.size() == mask2.size());
              for (std::size_t j = 0; j < m.size(); ++j) {
                m[j] = mask1[j] && mask2[j];
              }

              // Do the profile fitting
              ProfileFitter<double> fit(c, b, m.const_ref(), p, 1e-3, 100);
              // DIALS_ASSERT(fit.niter() < 100);

              // Set the data in the reflection
              intensity_val[i] = fit.intensity()[0];
              intensity_var[i] = fit.variance()[0];
              reference_cor[i] = fit.correlation();
              // reference_rmsd[i] = fit.rmsd();

              // Set the integrated flag
              flags[i] |= af::IntegratedPrf;
              success[i] = true;
            }

          } catch (dials::error const& e) {
            /* std::cout << e.what() << std::endl; */
            continue;
          }
        }
      }

      return success;
    }

    /**
     * Return a profile fitter
     * @return The profile fitter class
     */
    af::shared<bool> fit_detector_space(af::reflection_table reflections) const {
      // Check input is OK
      DIALS_ASSERT(reflections.is_consistent());
      DIALS_ASSERT(reflections.contains("shoebox"));
      DIALS_ASSERT(reflections.contains("flags"));
      DIALS_ASSERT(reflections.contains("partiality"));
      DIALS_ASSERT(reflections.contains("s1"));
      DIALS_ASSERT(reflections.contains("xyzcal.px"));
      DIALS_ASSERT(reflections.contains("xyzcal.mm"));

      // Get some data
      af::const_ref<Shoebox<>> sbox = reflections["shoebox"];
      af::const_ref<vec3<double>> s1 = reflections["s1"];
      af::const_ref<vec3<double>> xyzpx = reflections["xyzcal.px"];
      af::const_ref<vec3<double>> xyzmm = reflections["xyzcal.mm"];
      af::ref<std::size_t> flags = reflections["flags"];
      af::ref<double> intensity_val = reflections["intensity.prf.value"];
      af::ref<double> intensity_var = reflections["intensity.prf.variance"];
      af::ref<double> reference_cor = reflections["profile.correlation"];

      // Loop through all the reflections and process them
      af::shared<bool> success(reflections.size(), false);
      for (std::size_t i = 0; i < reflections.size(); ++i) {
        DIALS_ASSERT(sbox[i].is_consistent());

        // Set values to bad
        intensity_val[i] = 0.0;
        intensity_var[i] = -1.0;
        reference_cor[i] = 0.0;
        flags[i] &= ~af::IntegratedPrf;

        // Check if we want to use this reflection
        if (check2(flags[i], sbox[i])) {
          try {
            // Get the reference profiles
            std::size_t index = sampler_->nearest(sbox[i].panel, xyzpx[i]);
            data_const_reference d = data(index).const_ref();

            // Create the coordinate system
            vec3<double> m2 = spec_.goniometer().get_rotation_axis();
            vec3<double> s0 = spec_.beam()->get_s0();
            CoordinateSystem cs(m2, s0, s1[i], xyzmm[i][2]);

            // Compute the transform
            TransformReverse transform(spec_, cs, sbox[i].bbox, sbox[i].panel, d);

            // Get the transformed shoebox
            data_const_reference p = transform.profile().const_ref();

            // Create the data array
            af::versa<double, af::c_grid<3>> c(sbox[i].data.accessor());
            std::copy(sbox[i].data.begin(), sbox[i].data.end(), c.begin());

            // Create the background array
            af::versa<double, af::c_grid<3>> b(sbox[i].background.accessor());
            std::copy(sbox[i].background.begin(), sbox[i].background.end(), b.begin());

            // Create the mask array
            af::versa<bool, af::c_grid<3>> m(sbox[i].mask.accessor());

            std::transform(sbox[i].mask.begin(),
                           sbox[i].mask.end(),
                           m.begin(),
                           detail::check_mask_code(Valid | Foreground));

            // Do the profile fitting
            ProfileFitter<double> fit(
              c.const_ref(), b.const_ref(), m.const_ref(), p, 1e-3, 100);
            // DIALS_ASSERT(fit.niter() < 100);

            // Set the data in the reflection
            intensity_val[i] = fit.intensity()[0];
            intensity_var[i] = fit.variance()[0];
            reference_cor[i] = fit.correlation();

            // Set the integrated flag
            flags[i] |= af::IntegratedPrf;
            success[i] = true;

          } catch (dials::error const& e) {
            continue;
          }
        }
      }
      return success;
    }

    /**
     * @return a copy of the profile modller
     */
    pointer copy() const {
      GaussianRSProfileModeller result(beam_,
                                       detector_,
                                       goniometer_,
                                       scan_,
                                       sigma_b_,
                                       sigma_m_,
                                       n_sigma_,
                                       grid_size_,
                                       num_scan_points_,
                                       threshold_,
                                       grid_method_,
                                       fit_method_);
      result.finalized_ = finalized_;
      result.n_reflections_.assign(n_reflections_.begin(), n_reflections_.end());
      for (std::size_t i = 0; i < data_.size(); ++i) {
        if (data_[i].size() > 0) {
          result.data_[i] = data_type(accessor_, 0);
          result.mask_[i] = mask_type(accessor_, true);
          std::copy(data_[i].begin(), data_[i].end(), result.data_[i].begin());
          std::copy(mask_[i].begin(), mask_[i].end(), result.mask_[i].begin());
        }
      }
      return pointer(new GaussianRSProfileModeller(result));
    }

    void normalize_profiles() {
      finalize();
    }

  private:
    /**
     * Do we want to use the reflection in profile modelling
     * @param flags The reflection flags
     * @param partiality The reflection partiality
     * @param sbox The reflection shoebox
     * @return True/False
     */
    bool check1(std::size_t flags, double partiality, const Shoebox<>& sbox) const {
      // Check we're fully recorded
      bool full = partiality > 0.99;

      // Check reflection has been integrated
      bool integrated = flags & af::IntegratedSum;

      // Check if the bounding box is in the image
      bool bbox_valid = check_bbox_valid(flags, sbox);

      // Check if all pixels are valid
      bool pixels_valid = check_foreground_valid(flags, sbox);

      // Return whether to use or not
      return full && integrated && bbox_valid && pixels_valid;
    }

    /**
     * Do we want to use the reflection in profile fitting
     * @param flags The reflection flags
     * @param sbox The reflection shoebox
     * @return True/False
     */
    bool check2(std::size_t flags, const Shoebox<>& sbox) const {
      // Check if we want to integrate
      bool integrate = !(flags & af::DontIntegrate);

      // Check if the bounding box is in the image
      bool bbox_valid = check_bbox_valid(flags, sbox);

      // Check if all pixels are valid
      bool pixels_valid = check_foreground_valid(flags, sbox);

      // Return whether to use or not
      return integrate && bbox_valid && pixels_valid;
    }

    /**
     * Do we want to use the reflection in profile fitting
     * @param flags The reflection flags
     * @param sbox The reflection shoebox
     * @return True/False
     */
    bool check3(std::size_t flags, const Shoebox<>& sbox) const {
      // Check if we want to integrate
      bool integrate = !(flags & af::DontIntegrate);

      // Check if the bounding box is in the image
      bool bbox_valid = check_bbox_valid(flags, sbox);

      // Return whether to use or not
      return integrate && bbox_valid;
    }

    /**
     * Check if the bounding box is in entirely within the image
     * @param flags The reflection flags
     * @param sbox The reflection shoebox
     * @return True/False
     */
    bool check_bbox_valid(std::size_t flags, const Shoebox<>& sbox) const {
      return sbox.bbox[0] >= 0 && sbox.bbox[2] >= 0
             && sbox.bbox[1] <= spec_.detector()[sbox.panel].get_image_size()[0]
             && sbox.bbox[3] <= spec_.detector()[sbox.panel].get_image_size()[1];
    }

    /**
     * Check if all foreground pixels are valid
     * @param flags The reflection flags
     * @param sbox The reflection shoebox
     * @return True/False
     */
    bool check_foreground_valid(std::size_t flags, const Shoebox<>& sbox) const {
      bool pixels_valid = true;
      for (std::size_t i = 0; i < sbox.mask.size(); ++i) {
        if (sbox.mask[i] & Foreground && !(sbox.mask[i] & Valid)) {
          pixels_valid = false;
          break;
        }
      }
      return pixels_valid;
    }

    /**
     * Pre-evaluated reference profile in detector space for one sampler cell.
     */
    struct CellCube {
      af::versa<double, af::c_grid<3>> cube;  // (D_ext, H_ext, W_ext)
      af::versa<bool, af::c_grid<3>> mask;    // same shape
      scitbx::af::int6 bbox_nom;              // nominal bbox used to build the cube
      vec3<double> s1_nom;                    // nominal s1 for this cell
      double phi_nom;                         // nominal phi for this cell
      std::size_t panel;                      // panel index
    };

    /**
     * Result of extracting a shoebox-shaped sub-volume from a CellCube.
     */
    struct ExtractedRef {
      af::versa<double, af::c_grid<3>> profile;  // shoebox shape (D, H, W)
      af::versa<bool, af::c_grid<3>> mask;       // shoebox shape (D, H, W)
    };

    /**
     * Build a detector-space reference cube for the given sampler cell.
     *
     * For each voxel in the extended cube, computes the corresponding
     * reciprocal-space grid position and trilinearly interpolates the
     * 11x11x11 reference profile.
     */
    CellCube build_cell_cube(std::size_t cell_id, std::size_t panel) const {
      // Get spec parameters
      const vec3<double> s0 = spec_.beam()->get_s0();
      const vec3<double> m2 = spec_.goniometer().get_rotation_axis();
      const Scan& scan = spec_.scan();
      const auto& det_panel = spec_.detector()[panel];
      const double3 step_size = spec_.step_size();
      const double3 grid_cent = spec_.grid_centre();
      const int gs = static_cast<int>(spec_.half_grid_size());
      const int ref_size = 2 * gs + 1;

      // 1. Nominal centroid from sampler
      vec3<double> cell_coord = sampler_->coord(cell_id);
      double x_px = cell_coord[0];
      double y_px = cell_coord[1];
      double z_frame = cell_coord[2];

      // 2. Nominal s1: lab coord at pixel position, normalized to |s0|
      vec3<double> s1_nom = det_panel.get_pixel_lab_coord(vec2<double>(x_px, y_px));
      s1_nom = s1_nom.normalize() * s0.length();

      // 3. Nominal phi from scan
      double phi_nom = scan.get_angle_from_array_index(z_frame);

      // 4. Coordinate system at nominal position
      CoordinateSystem cs(m2, s0, s1_nom, phi_nom);

      // 5. Compute nominal bbox (same approach as BBoxCalculator3D)
      double delta_d = spec_.sigma_b() * spec_.n_sigma();
      double delta_m = spec_.sigma_m() * spec_.n_sigma();

      vec3<double> sd1 = cs.to_beam_vector(vec2<double>(-delta_d, -delta_d));
      vec3<double> sd2 = cs.to_beam_vector(vec2<double>(+delta_d, -delta_d));
      vec3<double> sd3 = cs.to_beam_vector(vec2<double>(-delta_d, +delta_d));
      vec3<double> sd4 = cs.to_beam_vector(vec2<double>(+delta_d, +delta_d));

      vec2<double> xy1 = det_panel.get_ray_intersection_px(sd1);
      vec2<double> xy2 = det_panel.get_ray_intersection_px(sd2);
      vec2<double> xy3 = det_panel.get_ray_intersection_px(sd3);
      vec2<double> xy4 = det_panel.get_ray_intersection_px(sd4);

      double phi1 = cs.to_rotation_angle_fast(-delta_m);
      double phi2 = cs.to_rotation_angle_fast(+delta_m);
      double z1 = scan.get_array_index_from_angle(phi1);
      double z2 = scan.get_array_index_from_angle(phi2);

      int x0 = static_cast<int>(
        std::floor(std::min(std::min(xy1[0], xy2[0]), std::min(xy3[0], xy4[0]))));
      int x1 = static_cast<int>(
        std::ceil(std::max(std::max(xy1[0], xy2[0]), std::max(xy3[0], xy4[0]))));
      int y0 = static_cast<int>(
        std::floor(std::min(std::min(xy1[1], xy2[1]), std::min(xy3[1], xy4[1]))));
      int y1 = static_cast<int>(
        std::ceil(std::max(std::max(xy1[1], xy2[1]), std::max(xy3[1], xy4[1]))));
      int z0 = static_cast<int>(std::floor(std::min(z1, z2)));
      int z1i = static_cast<int>(std::ceil(std::max(z1, z2)));

      scitbx::af::int6 bbox_nom;
      bbox_nom[0] = x0;
      bbox_nom[1] = x1;
      bbox_nom[2] = y0;
      bbox_nom[3] = y1;
      bbox_nom[4] = z0;
      bbox_nom[5] = z1i;

      // 6. Extended cube dimensions (+1 on each side for trilinear margin)
      int W_sbox = x1 - x0;
      int H_sbox = y1 - y0;
      int D_sbox = z1i - z0;
      int W_ext = W_sbox + 2;
      int H_ext = H_sbox + 2;
      int D_ext = D_sbox + 2;

      af::c_grid<3> cube_acc(D_ext, H_ext, W_ext);
      af::versa<double, af::c_grid<3>> cube(cube_acc, 0.0);
      af::versa<bool, af::c_grid<3>> cube_mask(cube_acc, false);

      // Reference profile for this cell: layout (gj, gi, gk)
      data_const_reference ref = data(cell_id).const_ref();
      mask_const_reference ref_mask = mask(cell_id).const_ref();

      // Precompute e1/e2 scaled by 1/|s1| (same as TransformForward)
      vec3<double> e1_scaled = cs.e1_axis() / s1_nom.length();
      vec3<double> e2_scaled = cs.e2_axis() / s1_nom.length();

      // 7. For each voxel in the extended cube
      for (int d = 0; d < D_ext; ++d) {
        // Frame index (the -1 accounts for the +1 extension margin)
        double z_vox = static_cast<double>(z0 + d - 1);
        double phi_vox = scan.get_angle_from_array_index(z_vox);
        double c3 = cs.from_rotation_angle_fast(phi_vox);
        double gk = grid_cent[0] + c3 / step_size[0];

        for (int h = 0; h < H_ext; ++h) {
          double y_vox = static_cast<double>(y0 + h - 1);

          for (int w = 0; w < W_ext; ++w) {
            double x_vox = static_cast<double>(x0 + w - 1);

            // Beam vector at this pixel
            vec3<double> s1_vox =
              det_panel.get_pixel_lab_coord(vec2<double>(x_vox, y_vox));
            s1_vox = s1_vox.normalize() * s0.length();

            // Kabsch e1, e2 coords
            vec3<double> ds = s1_vox - s1_nom;
            double gi = grid_cent[2] + (e1_scaled * ds) / step_size[2];
            double gj = grid_cent[1] + (e2_scaled * ds) / step_size[1];

            // Trilinear interpolation of the reference profile
            // Integer base indices
            int ik = static_cast<int>(std::floor(gk));
            int ij = static_cast<int>(std::floor(gj));
            int ii = static_cast<int>(std::floor(gi));

            // Check all 8 neighbors are in bounds
            if (ik < 0 || ik + 1 >= ref_size || ij < 0 || ij + 1 >= ref_size || ii < 0
                || ii + 1 >= ref_size) {
              cube(d, h, w) = 0.0;
              cube_mask(d, h, w) = false;
              continue;
            }

            // Check reference mask at all 8 corners
            bool all_mask = ref_mask(ik, ij, ii) && ref_mask(ik + 1, ij, ii)
                            && ref_mask(ik, ij, ii + 1) && ref_mask(ik + 1, ij, ii + 1)
                            && ref_mask(ik, ij + 1, ii) && ref_mask(ik + 1, ij + 1, ii)
                            && ref_mask(ik, ij + 1, ii + 1)
                            && ref_mask(ik + 1, ij + 1, ii + 1);

            if (!all_mask) {
              cube(d, h, w) = 0.0;
              cube_mask(d, h, w) = false;
              continue;
            }

            // Fractional parts
            double fk = gk - ik;
            double fj = gj - ij;
            double fi = gi - ii;

            // Trilinear interpolation: ref is indexed as (e3, e2, e1) = (gk, gj, gi)
            double val = ref(ik, ij, ii) * (1 - fi) * (1 - fj) * (1 - fk)
                         + ref(ik + 1, ij, ii) * (1 - fi) * (1 - fj) * fk
                         + ref(ik, ij, ii + 1) * fi * (1 - fj) * (1 - fk)
                         + ref(ik + 1, ij, ii + 1) * fi * (1 - fj) * fk
                         + ref(ik, ij + 1, ii) * (1 - fi) * fj * (1 - fk)
                         + ref(ik + 1, ij + 1, ii) * (1 - fi) * fj * fk
                         + ref(ik, ij + 1, ii + 1) * fi * fj * (1 - fk)
                         + ref(ik + 1, ij + 1, ii + 1) * fi * fj * fk;

            cube(d, h, w) = val;
            cube_mask(d, h, w) = true;
          }
        }
      }

      CellCube result;
      result.cube = cube;
      result.mask = cube_mask;
      result.bbox_nom = bbox_nom;
      result.s1_nom = s1_nom;
      result.phi_nom = phi_nom;
      result.panel = panel;
      return result;
    }

    /**
     * Extract a shoebox-shaped sub-volume from a pre-computed CellCube via
     * trilinear interpolation at a given sub-pixel offset.
     *
     * @param cc   The pre-computed cell cube (extended by +1 on each side)
     * @param dx   Sub-pixel offset in x (columns, pixels)
     * @param dy   Sub-pixel offset in y (rows, pixels)
     * @param dz   Sub-pixel offset in z (frames)
     * @param bbox The reflection bounding box [x0,x1,y0,y1,z0,z1]
     * @return ExtractedRef with profile and mask of shoebox shape
     */
    ExtractedRef extract_from_cube(const CellCube& cc,
                                   double dx,
                                   double dy,
                                   double dz,
                                   int6 bbox) const {
      // Output shoebox dimensions from bbox
      int W = bbox[1] - bbox[0];
      int H = bbox[3] - bbox[2];
      int D = bbox[5] - bbox[4];

      // Cube extended dimensions
      int W_ext = static_cast<int>(cc.cube.accessor()[2]);
      int H_ext = static_cast<int>(cc.cube.accessor()[1]);
      int D_ext = static_cast<int>(cc.cube.accessor()[0]);

      // Allocate output arrays
      af::c_grid<3> out_acc(D, H, W);
      af::versa<double, af::c_grid<3>> result_profile(out_acc, 0.0);
      af::versa<bool, af::c_grid<3>> result_mask(out_acc, false);

      for (int d = 0; d < D; ++d) {
        for (int h = 0; h < H; ++h) {
          for (int w = 0; w < W; ++w) {
            // Position in the extended cube (+1.0 enters the margin zone)
            double cx = w + 1.0 + dx;
            double cy = h + 1.0 + dy;
            double cz = d + 1.0 + dz;

            // Integer base and fractional parts
            int ix = static_cast<int>(std::floor(cx));
            int iy = static_cast<int>(std::floor(cy));
            int iz = static_cast<int>(std::floor(cz));
            double fx = cx - ix;
            double fy = cy - iy;
            double fz = cz - iz;

            // Bounds check: all 8 corners must be in [0, dim-1]
            bool in_bounds = ix >= 0 && ix + 1 < W_ext && iy >= 0 && iy + 1 < H_ext
                             && iz >= 0 && iz + 1 < D_ext;

            if (!in_bounds) {
              result_profile(d, h, w) = 0.0;
              result_mask(d, h, w) = false;
              continue;
            }

            // Trilinear interpolation (8 loads, 7 muls, 7 adds)
            double val = cc.cube(iz, iy, ix) * (1 - fx) * (1 - fy) * (1 - fz)
                         + cc.cube(iz, iy, ix + 1) * fx * (1 - fy) * (1 - fz)
                         + cc.cube(iz, iy + 1, ix) * (1 - fx) * fy * (1 - fz)
                         + cc.cube(iz, iy + 1, ix + 1) * fx * fy * (1 - fz)
                         + cc.cube(iz + 1, iy, ix) * (1 - fx) * (1 - fy) * fz
                         + cc.cube(iz + 1, iy, ix + 1) * fx * (1 - fy) * fz
                         + cc.cube(iz + 1, iy + 1, ix) * (1 - fx) * fy * fz
                         + cc.cube(iz + 1, iy + 1, ix + 1) * fx * fy * fz;

            result_profile(d, h, w) = val;

            // Mask: true only if ALL 8 cube mask corners are true
            result_mask(d, h, w) =
              cc.mask(iz, iy, ix) && cc.mask(iz, iy, ix + 1) && cc.mask(iz, iy + 1, ix)
              && cc.mask(iz, iy + 1, ix + 1) && cc.mask(iz + 1, iy, ix)
              && cc.mask(iz + 1, iy, ix + 1) && cc.mask(iz + 1, iy + 1, ix)
              && cc.mask(iz + 1, iy + 1, ix + 1);
          }
        }
      }

      ExtractedRef extracted;
      extracted.profile = result_profile;
      extracted.mask = result_mask;
      return extracted;
    }

    TransformSpec spec_;

    // Cell-cache: per-cell precomputed detector-space reference cubes
    mutable std::vector<CellCube> cube_cache_;
    mutable std::vector<bool> cube_valid_;
  };

}}  // namespace dials::algorithms

#endif  // DIALS_ALGORITHMS_PROFILE_MODEL_GAUSSIAN_RS_MODELLER_H
