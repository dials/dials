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

#include <memory>
#include <fstream>
#include <limits>
#include <algorithm>
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

    enum FitMethod { ReciprocalSpace = 1, DetectorSpace = 2 };

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
      uint8_t mask_code;
      check_mask_code(uint8_t code) : mask_code(code) {}
      bool operator()(uint8_t a) const {
        return ((a & mask_code) == mask_code);
      }
    };

    struct check_either_mask_code {
      uint8_t mask_code1;
      uint8_t mask_code2;
      check_either_mask_code(uint8_t code1, uint8_t code2)
          : mask_code1(code1), mask_code2(code2) {}
      bool operator()(uint8_t a) const {
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
      af::const_ref<Shoebox<> > sbox = reflections["shoebox"];
      af::const_ref<double> partiality = reflections["partiality"];
      af::const_ref<vec3<double> > s1 = reflections["s1"];
      af::const_ref<vec3<double> > xyzpx = reflections["xyzcal.px"];
      af::const_ref<vec3<double> > xyzmm = reflections["xyzcal.mm"];
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
          af::versa<double, af::c_grid<3> > data(sbox[i].data.accessor());
          std::transform(sbox[i].data.begin(),
                         sbox[i].data.end(),
                         sbox[i].background.begin(),
                         data.begin(),
                         std::minus<double>());

          // Create the mask array
          af::versa<bool, af::c_grid<3> > mask(sbox[i].mask.accessor());
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
     * A reflection's shoebox carried onto the grid of the reference profiles,
     * ready to be fitted. The transform does not depend on the reference
     * profiles, so it can be computed as soon as the shoebox is complete and the
     * fit done once the profile it needs has been finalized.
     */
    struct PreparedFit {
      std::size_t index;  // the reference profile to fit against
      data_type profile;
      data_type background;
      mask_type mask;
    };

    struct FitResult {
      double intensity;
      double variance;
      double correlation;
    };

    /**
     * Transform a shoebox onto the reference profile grid for fitting in
     * reciprocal space: the first half of fit_reciprocal_space.
     */
    PreparedFit prepare_fit_reciprocal_space(const Shoebox<>& sbox,
                                             const vec3<double>& s1,
                                             const vec3<double>& xyzpx,
                                             const vec3<double>& xyzmm) const {
      PreparedFit result;

      // Get the reference profile index
      result.index = sampler_->nearest(sbox.panel, xyzpx);

      // Create the coordinate system
      vec3<double> m2 = spec_.goniometer().get_rotation_axis();
      vec3<double> s0 = spec_.beam()->get_s0();
      CoordinateSystem cs(m2, s0, s1, xyzmm[2]);

      // Create the data array
      af::versa<double, af::c_grid<3> > data(sbox.data.accessor());
      std::copy(sbox.data.begin(), sbox.data.end(), data.begin());

      // Create the background array
      af::versa<double, af::c_grid<3> > background(sbox.background.accessor());
      std::copy(sbox.background.begin(), sbox.background.end(), background.begin());

      // Create the mask array
      af::versa<bool, af::c_grid<3> > mask(sbox.mask.accessor());
      std::transform(sbox.mask.begin(),
                     sbox.mask.end(),
                     mask.begin(),
                     detail::check_mask_code(Valid | Foreground));

      // Compute the transform
      TransformForward<double> transform(spec_,
                                         cs,
                                         sbox.bbox,
                                         sbox.panel,
                                         data.const_ref(),
                                         background.const_ref(),
                                         mask.const_ref());

      // Keep the transformed shoebox
      result.profile = transform.profile();
      result.background = transform.background();
      result.mask = transform.mask();
      return result;
    }

    /**
     * Fit a transformed shoebox against its reference profile: the second half
     * of fit_reciprocal_space. Throws if the profile is empty.
     */
    FitResult fit_prepared_reciprocal_space(const PreparedFit& prepared) const {
      // Get the reference profile
      data_const_reference p = data(prepared.index).const_ref();
      mask_const_reference mask1 = mask(prepared.index).const_ref();

      // Combine the masks
      mask_const_reference mask2 = prepared.mask.const_ref();
      af::versa<bool, af::c_grid<3> > m(mask2.accessor());
      DIALS_ASSERT(mask1.size() == mask2.size());
      for (std::size_t j = 0; j < m.size(); ++j) {
        m[j] = mask1[j] && mask2[j];
      }

      // Do the profile fitting
      ProfileFitter<double> fit(prepared.profile.const_ref(),
                                prepared.background.const_ref(),
                                m.const_ref(),
                                p,
                                1e-3,
                                100);
      // DIALS_ASSERT(fit.niter() < 100);

      FitResult result;
      result.intensity = fit.intensity()[0];
      result.variance = fit.variance()[0];
      result.correlation = fit.correlation();
      return result;
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
      af::const_ref<Shoebox<> > sbox = reflections["shoebox"];
      af::const_ref<vec3<double> > s1 = reflections["s1"];
      af::const_ref<vec3<double> > xyzpx = reflections["xyzcal.px"];
      af::const_ref<vec3<double> > xyzmm = reflections["xyzcal.mm"];
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
            // Transform the shoebox and fit it against its reference profile
            PreparedFit prepared =
              prepare_fit_reciprocal_space(sbox[i], s1[i], xyzpx[i], xyzmm[i]);
            FitResult fit = fit_prepared_reciprocal_space(prepared);

            // Set the data in the reflection
            intensity_val[i] = fit.intensity;
            intensity_var[i] = fit.variance;
            reference_cor[i] = fit.correlation;
            // reference_rmsd[i] = fit.rmsd();

            // Set the integrated flag
            flags[i] |= af::IntegratedPrf;
            success[i] = true;

          } catch (dials::error const& e) {
            /* std::cout << e.what() << std::endl; */
            continue;
          }
        }
      }
      return success;
    }

    /**
     * For each reference profile, the last frame (exclusive) of any of these
     * reflections that could contribute to it in model(): after that frame has
     * been processed, nothing more can be added to the profile. The minimum
     * int where none can.
     */
    af::shared<int> learning_deadlines(af::reflection_table reflections) const {
      DIALS_ASSERT(reflections.contains("bbox"));
      DIALS_ASSERT(reflections.contains("panel"));
      DIALS_ASSERT(reflections.contains("xyzcal.px"));
      af::const_ref<int6> bbox = reflections["bbox"];
      af::const_ref<std::size_t> panel = reflections["panel"];
      af::const_ref<vec3<double> > xyzpx = reflections["xyzcal.px"];
      af::shared<int> result(sampler_->size(), std::numeric_limits<int>::min());
      for (std::size_t i = 0; i < reflections.size(); ++i) {
        af::shared<std::size_t> indices;
        try {
          indices = sampler_->nearest_n(panel[i], xyzpx[i]);
        } catch (dials::error const&) {
          continue;
        }
        for (std::size_t j = 0; j < indices.size(); ++j) {
          DIALS_ASSERT(indices[j] < result.size());
          result[indices[j]] = std::max(result[indices[j]], bbox[i][5]);
        }
      }
      return result;
    }

    /**
     * For each reflection, the reference profile it is fitted against, or -1
     * where there is none.
     */
    af::shared<int> fitting_cells(af::reflection_table reflections) const {
      DIALS_ASSERT(reflections.contains("panel"));
      DIALS_ASSERT(reflections.contains("xyzcal.px"));
      af::const_ref<std::size_t> panel = reflections["panel"];
      af::const_ref<vec3<double> > xyzpx = reflections["xyzcal.px"];
      af::shared<int> result(reflections.size(), -1);
      for (std::size_t i = 0; i < reflections.size(); ++i) {
        try {
          result[i] = (int)sampler_->nearest(panel[i], xyzpx[i]);
        } catch (dials::error const&) {
          continue;
        }
      }
      return result;
    }

    /**
     * For each reflection, whether model() could add it to any of the selected
     * reference profiles.
     * @param reflections The reflections
     * @param cells A selection of the reference profiles
     */
    af::shared<bool> contributes_to(af::reflection_table reflections,
                                    const af::const_ref<bool>& cells) const {
      DIALS_ASSERT(reflections.contains("panel"));
      DIALS_ASSERT(reflections.contains("xyzcal.px"));
      DIALS_ASSERT(cells.size() == sampler_->size());
      af::const_ref<std::size_t> panel = reflections["panel"];
      af::const_ref<vec3<double> > xyzpx = reflections["xyzcal.px"];
      af::shared<bool> result(reflections.size(), false);
      for (std::size_t i = 0; i < reflections.size(); ++i) {
        af::shared<std::size_t> indices;
        try {
          indices = sampler_->nearest_n(panel[i], xyzpx[i]);
        } catch (dials::error const&) {
          continue;
        }
        for (std::size_t j = 0; j < indices.size(); ++j) {
          DIALS_ASSERT(indices[j] < cells.size());
          if (cells[indices[j]]) {
            result[i] = true;
            break;
          }
        }
      }
      return result;
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
      af::const_ref<Shoebox<> > sbox = reflections["shoebox"];
      af::const_ref<vec3<double> > s1 = reflections["s1"];
      af::const_ref<vec3<double> > xyzpx = reflections["xyzcal.px"];
      af::const_ref<vec3<double> > xyzmm = reflections["xyzcal.mm"];
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
            af::versa<double, af::c_grid<3> > c(sbox[i].data.accessor());
            std::copy(sbox[i].data.begin(), sbox[i].data.end(), c.begin());

            // Create the background array
            af::versa<double, af::c_grid<3> > b(sbox[i].background.accessor());
            std::copy(sbox[i].background.begin(), sbox[i].background.end(), b.begin());

            // Create the mask array
            af::versa<bool, af::c_grid<3> > m(sbox[i].mask.accessor());

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
      result.cell_finalized_ = cell_finalized_;
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

    TransformSpec spec_;
  };

  /**
   * Profile fits waiting for their reference profiles, for integrating in one
   * pass over the images. A reflection's shoebox is transformed onto the
   * reference profile grid as soon as it is complete and only the transform is
   * kept; it is fitted, and the transform released, once its reference profile
   * has been finalized. Every outcome is recorded by row, success or failure,
   * exactly as fit_reciprocal_space would have given it.
   */
  class PendingFits {
  public:
    typedef GaussianRSProfileModeller::PreparedFit PreparedFit;
    typedef GaussianRSProfileModeller::FitResult FitResult;

    PendingFits() : held_bytes_(0), peak_held_(0), peak_bytes_(0) {}

    /**
     * Prepare the fits of some reflections and hold them
     * @param modeller The profile modeller for their experiment
     * @param reflections The reflections, with complete shoeboxes
     * @param rows An identifier for each reflection, used in the results
     */
    void add(const GaussianRSProfileModeller& modeller,
             af::reflection_table reflections,
             const af::const_ref<std::size_t>& rows) {
      DIALS_ASSERT(reflections.is_consistent());
      DIALS_ASSERT(reflections.contains("shoebox"));
      DIALS_ASSERT(reflections.contains("flags"));
      DIALS_ASSERT(reflections.contains("s1"));
      DIALS_ASSERT(reflections.contains("xyzcal.px"));
      DIALS_ASSERT(reflections.contains("xyzcal.mm"));
      DIALS_ASSERT(rows.size() == reflections.size());
      af::const_ref<Shoebox<> > sbox = reflections["shoebox"];
      af::const_ref<vec3<double> > s1 = reflections["s1"];
      af::const_ref<vec3<double> > xyzpx = reflections["xyzcal.px"];
      af::const_ref<vec3<double> > xyzmm = reflections["xyzcal.mm"];
      af::const_ref<std::size_t> flags = reflections["flags"];
      for (std::size_t i = 0; i < reflections.size(); ++i) {
        DIALS_ASSERT(sbox[i].is_consistent());
        if (flags[i] & af::DontIntegrate) {
          record_failure(rows[i]);
          continue;
        }
        try {
          Pending pending;
          pending.row = rows[i];
          pending.fit =
            modeller.prepare_fit_reciprocal_space(sbox[i], s1[i], xyzpx[i], xyzmm[i]);
          held_bytes_ += bytes(pending.fit);
          pending_.push_back(pending);
        } catch (dials::error const&) {
          record_failure(rows[i]);
        }
      }
      peak_held_ = std::max(peak_held_, pending_.size());
      peak_bytes_ = std::max(peak_bytes_, held_bytes_);
    }

    /**
     * Fit every held reflection whose reference profile is finalized, and
     * release it
     * @param modeller The profile modeller the reflections were prepared with
     * @return The number fitted
     */
    std::size_t fit_ready(const GaussianRSProfileModeller& modeller) {
      std::size_t kept = 0;
      std::size_t fitted = 0;
      for (std::size_t i = 0; i < pending_.size(); ++i) {
        if (modeller.cell_finalized(pending_[i].fit.index)) {
          try {
            FitResult r = modeller.fit_prepared_reciprocal_space(pending_[i].fit);
            rows_.push_back(pending_[i].row);
            intensity_.push_back(r.intensity);
            variance_.push_back(r.variance);
            correlation_.push_back(r.correlation);
            success_.push_back(true);
          } catch (dials::error const&) {
            record_failure(pending_[i].row);
          }
          held_bytes_ -= bytes(pending_[i].fit);
          fitted++;
        } else {
          if (kept != i) {
            pending_[kept] = pending_[i];
          }
          kept++;
        }
      }
      pending_.resize(kept);
      return fitted;
    }

    std::size_t held() const {
      return pending_.size();
    }

    std::size_t held_bytes() const {
      return held_bytes_;
    }

    std::size_t peak_held() const {
      return peak_held_;
    }

    std::size_t peak_bytes() const {
      return peak_bytes_;
    }

    af::shared<std::size_t> rows() const {
      return rows_;
    }

    af::shared<double> intensity() const {
      return intensity_;
    }

    af::shared<double> variance() const {
      return variance_;
    }

    af::shared<double> correlation() const {
      return correlation_;
    }

    af::shared<bool> success() const {
      return success_;
    }

  private:
    struct Pending {
      std::size_t row;
      PreparedFit fit;
    };

    // The values fit_reciprocal_space leaves on a reflection it cannot fit
    void record_failure(std::size_t row) {
      rows_.push_back(row);
      intensity_.push_back(0.0);
      variance_.push_back(-1.0);
      correlation_.push_back(0.0);
      success_.push_back(false);
    }

    static std::size_t bytes(const PreparedFit& fit) {
      return (fit.profile.size() + fit.background.size()) * sizeof(double)
             + fit.mask.size() * sizeof(bool);
    }

    std::vector<Pending> pending_;
    std::size_t held_bytes_;
    std::size_t peak_held_;
    std::size_t peak_bytes_;
    af::shared<std::size_t> rows_;
    af::shared<double> intensity_;
    af::shared<double> variance_;
    af::shared<double> correlation_;
    af::shared<bool> success_;
  };

}}  // namespace dials::algorithms

#endif  // DIALS_ALGORITHMS_PROFILE_MODEL_GAUSSIAN_RS_MODELLER_H
