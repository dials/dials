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
#include <vector>
#include <dials/algorithms/profile_model/gaussian_rs/transform/transform.h>
#include <dials/algorithms/profile_model/modeller/empirical_modeller.h>
#include <dials/algorithms/profile_model/modeller/single_sampler.h>
#include <dials/algorithms/profile_model/modeller/grid_sampler.h>
#include <dials/algorithms/profile_model/modeller/circle_sampler.h>
#include <dials/algorithms/profile_model/modeller/ewald_sphere_sampler.h>
#include <dials/algorithms/integration/fit/fitting.h>

namespace dials { namespace algorithms {

  using dials::algorithms::profile_model::gaussian_rs::CoordinateSystem;
  using dials::algorithms::profile_model::gaussian_rs::transform::MapFramesForward;
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
              // --- CELL-CACHE PATH: linearized gc() mapping into reference ---
              if (!valid(index)) continue;

              // Reference profile and mask (11x11x11, layout: (e2, e1, e3))
              data_const_reference ref = data(index).const_ref();
              mask_const_reference ref_mask = mask(index).const_ref();
              const int n_grid = static_cast<int>(2 * spec_.half_grid_size() + 1);

              // Construct CoordinateSystem at reflection position
              vec3<double> m2 = spec_.goniometer().get_rotation_axis();
              vec3<double> s0 = spec_.beam()->get_s0();
              CoordinateSystem cs(m2, s0, s1[i], xyzmm[i][2]);

              // Detector geometry
              const Panel& panel = spec_.detector()[sbox[i].panel];
              vec3<double> d_fast = panel.get_fast_axis();
              vec3<double> d_slow = panel.get_slow_axis();
              vec2<double> px = panel.get_pixel_size();

              // Lab position and beam direction at reflection centroid
              vec3<double> p =
                panel.get_pixel_lab_coord(vec2<double>(xyzpx[i][0], xyzpx[i][1]));
              double r = p.length();
              vec3<double> p_hat = p / r;

              // Kabsch basis vectors scaled by 1/|s1|
              double s1_len = cs.s1().length();
              vec3<double> e1_s = cs.e1_axis() / s1_len;
              vec3<double> e2_s = cs.e2_axis() / s1_len;

              // Project panel axes perpendicular to beam direction
              vec3<double> f_perp = d_fast - p_hat * (p_hat * d_fast);
              vec3<double> s_perp = d_slow - p_hat * (p_hat * d_slow);

              // Scale factor
              double scale = s1_len / r;

              // 2x2 spatial Jacobian: d(c1,c2)/d(x_px,y_px)
              double j11 = scale * px[0] * (e1_s * f_perp);
              double j12 = scale * px[1] * (e1_s * s_perp);
              double j21 = scale * px[0] * (e2_s * f_perp);
              double j22 = scale * px[1] * (e2_s * s_perp);

              // Convert to grid units: divide by step size
              // step_size() returns (step_e3, step_e2, step_e1)
              double g00 = j11 / spec_.step_size()[2];  // dc1_grid/dx
              double g01 = j12 / spec_.step_size()[2];  // dc1_grid/dy
              double g10 = j21 / spec_.step_size()[1];  // dc2_grid/dx
              double g11 = j22 / spec_.step_size()[1];  // dc2_grid/dy

              // e3 (frame) direction
              double osc_rad = spec_.scan().get_oscillation()[1];
              double g22 = cs.zeta() * osc_rad / spec_.step_size()[0];

              // Grid center: reflection centroid maps to (0,0,0) in its
              // own Kabsch frame, which maps to grid center by construction
              double gc0 = spec_.grid_centre()[2];  // e1 grid center
              double gc1 = spec_.grid_centre()[1];  // e2 grid center
              double gc2 = spec_.grid_centre()[0];  // e3 grid center

              // Shoebox dimensions
              int6 bbox = sbox[i].bbox;
              int W = bbox[1] - bbox[0];
              int H = bbox[3] - bbox[2];
              int D = bbox[5] - bbox[4];

              // --- Hybrid reverse interpolation: bilinear spatial +
              //     MapFramesForward z_fraction along e3 ---
              // Invert the 2x2 spatial block of the Jacobian G once per
              // reflection.  The e3 direction is handled separately via the
              // Kabsch-2010 z_fraction matrix, so we no longer need inv22.
              double det2 = g00 * g11 - g01 * g10;

              const double DET_EPS = 1e-6;
              if (std::abs(det2) < DET_EPS || std::abs(g22) < DET_EPS) {
                // Ill-conditioned (zeta ~ 0 or degenerate spatial block);
                // skip this reflection rather than falling through.
                continue;
              }

              double inv00 = g11 / det2;
              double inv01 = -g01 / det2;
              double inv10 = -g10 / det2;
              double inv11 = g00 / det2;

              // Allocate grid-space buffers (11x11x11, (jj, ii, kk) layout
              // matching `ref`).
              af::c_grid<3> grid_acc(n_grid, n_grid, n_grid);
              af::versa<double, af::c_grid<3>> data_grid(grid_acc, 0.0);
              af::versa<double, af::c_grid<3>> bg_grid(grid_acc, 0.0);
              af::versa<bool, af::c_grid<3>> combined_mask(grid_acc, false);

              // Build the per-frame z_fraction matrix via MapFramesForward.
              // zfraction(k, kk) is the fraction of frame-k's angular mass
              // that lands in e3 grid cell kk; summing over k gives the
              // Kabsch partiality for cell kk.  This replaces the expensive
              // quad_to_grid polygon clipping used by TransformForward and
              // provides the e3 weighting we'd otherwise get from a dense
              // sub-sampling loop.
              MapFramesForward<double> map_fwd(spec_.scan().get_array_range()[0],
                                               spec_.scan().get_oscillation()[0],
                                               spec_.scan().get_oscillation()[1],
                                               spec_.sigma_m(),
                                               spec_.n_sigma(),
                                               spec_.half_grid_size());
              af::versa<double, af::c_grid<2>> zfraction_arr =
                map_fwd(vec2<int>(bbox[4], bbox[5]), cs.phi(), cs.zeta());
              af::const_ref<double, af::c_grid<2>> zfraction =
                zfraction_arr.const_ref();

              // Spatial supersampling only: 2D density is |det2| (grid
              // cells per pixel^2).  The e3 axis is now handled analytically
              // by z_fraction, so there is no sub-sampling in that direction.
              double sp_density = std::abs(det2);
              int N_sp = std::max(1, (int)std::ceil(std::sqrt(1.0 / sp_density)));
              double inv_Nsp = 1.0 / N_sp;

              // Guard against pathological spatial shapes.
              const int MAX_SPATIAL_SAMPLES = 64;
              if (N_sp * N_sp > MAX_SPATIAL_SAMPLES) {
                continue;
              }

              const int need = Valid | Foreground;

              // Per-frame bilinear stack buffers (length D, one entry per
              // shoebox frame).  Re-used across (jj, ii) cells.
              std::vector<double> V_data(D, 0.0);
              std::vector<double> V_bg(D, 0.0);
              std::vector<bool> V_fg_ok(D, true);

              // Walk the 121 spatial (jj, ii) cells.  For each one we
              // compute a bilinear interpolation at the reverse-mapped
              // shoebox column (vx, vy) for every frame k in the shoebox,
              // producing the per-frame stack vectors.  Then we collapse
              // the stack into e3 grid cells via the z_fraction weights.
              for (int jj = 0; jj < n_grid; ++jj) {
                for (int ii = 0; ii < n_grid; ++ii) {
                  std::fill(V_data.begin(), V_data.end(), 0.0);
                  std::fill(V_bg.begin(), V_bg.end(), 0.0);
                  std::fill(V_fg_ok.begin(), V_fg_ok.end(), true);
                  int n_hit_spatial = 0;

                  // Spatial supersampling (2D only).
                  for (int sa = 0; sa < N_sp; ++sa) {
                    double dgj = (jj + (sa + 0.5) * inv_Nsp) - gc1;
                    for (int sb = 0; sb < N_sp; ++sb) {
                      double dgi = (ii + (sb + 0.5) * inv_Nsp) - gc0;

                      // Inverse spatial map: grid offset -> pixel offset.
                      double dx = inv00 * dgi + inv01 * dgj;
                      double dy = inv10 * dgi + inv11 * dgj;
                      double x = xyzpx[i][0] + dx;
                      double y = xyzpx[i][1] + dy;

                      // Shoebox voxel coordinates (continuous, with the
                      // -0.5 pixel-centre convention).
                      double vx = x - bbox[0] - 0.5;
                      double vy = y - bbox[2] - 0.5;

                      int iw = static_cast<int>(std::floor(vx));
                      int ih = static_cast<int>(std::floor(vy));
                      double fw = vx - iw;
                      double fh = vy - ih;

                      // Symmetric edge clamp (same as pure-trilinear path).
                      if (iw < 0) {
                        iw = 0;
                        fw = 0.0;
                      }
                      if (ih < 0) {
                        ih = 0;
                        fh = 0.0;
                      }
                      if (iw >= W - 1) {
                        iw = W - 2;
                        fw = 1.0;
                      }
                      if (ih >= H - 1) {
                        ih = H - 2;
                        fh = 1.0;
                      }

                      // Bilinear weights.
                      double w00 = (1.0 - fh) * (1.0 - fw);
                      double w01 = (1.0 - fh) * fw;
                      double w10 = fh * (1.0 - fw);
                      double w11 = fh * fw;

                      // Bilinear interp at (vx, vy) for every frame k.
                      for (int k = 0; k < D; ++k) {
                        V_data[k] +=
                          w00 * static_cast<double>(sbox[i].data(k, ih, iw))
                          + w01 * static_cast<double>(sbox[i].data(k, ih, iw + 1))
                          + w10 * static_cast<double>(sbox[i].data(k, ih + 1, iw))
                          + w11 * static_cast<double>(sbox[i].data(k, ih + 1, iw + 1));
                        V_bg[k] +=
                          w00 * static_cast<double>(sbox[i].background(k, ih, iw))
                          + w01 * static_cast<double>(sbox[i].background(k, ih, iw + 1))
                          + w10 * static_cast<double>(sbox[i].background(k, ih + 1, iw))
                          + w11
                              * static_cast<double>(
                                sbox[i].background(k, ih + 1, iw + 1));
                        bool fg_ok =
                          ((sbox[i].mask(k, ih, iw) & need) == need)
                          && ((sbox[i].mask(k, ih, iw + 1) & need) == need)
                          && ((sbox[i].mask(k, ih + 1, iw) & need) == need)
                          && ((sbox[i].mask(k, ih + 1, iw + 1) & need) == need);
                        V_fg_ok[k] = V_fg_ok[k] && fg_ok;
                      }
                      ++n_hit_spatial;
                    }
                  }

                  if (n_hit_spatial == 0) continue;

                  // Average over spatial supersamples.
                  double inv_nh = 1.0 / static_cast<double>(n_hit_spatial);
                  for (int k = 0; k < D; ++k) {
                    V_data[k] *= inv_nh;
                    V_bg[k] *= inv_nh;
                  }

                  // Spatial cell-volume scaling: 1/|det2| ONLY (NOT 1/det_G).
                  // z_fraction is already a normalised fraction-of-mass and
                  // therefore accounts for the e3 axis volume on its own.
                  double inv_det2 = 1.0 / std::abs(det2);

                  // Collapse the per-frame stack into e3 grid cells via
                  // the z_fraction weights.
                  for (int kk = 0; kk < n_grid; ++kk) {
                    double accum_d = 0.0;
                    double accum_b = 0.0;
                    for (int k = 0; k < D; ++k) {
                      double zf = zfraction(k, kk);
                      accum_d += V_data[k] * zf;
                      accum_b += V_bg[k] * zf;
                    }
                    data_grid(jj, ii, kk) = accum_d * inv_det2;
                    bg_grid(jj, ii, kk) = accum_b * inv_det2;

                    // Mask: cell is valid only if the reference mask is true
                    // here AND every frame with a non-trivial z_fraction
                    // contribution has a fully-foreground bilinear stencil.
                    bool ok = ref_mask(jj, ii, kk);
                    if (ok) {
                      for (int k = 0; k < D; ++k) {
                        if (zfraction(k, kk) > 1e-6 && !V_fg_ok[k]) {
                          ok = false;
                          break;
                        }
                      }
                    }
                    combined_mask(jj, ii, kk) = ok;
                  }
                }
              }

              // Partiality guard: require a minimum number of valid cells.
              std::size_t num_valid = 0;
              for (std::size_t t = 0; t < combined_mask.size(); ++t) {
                if (combined_mask[t]) ++num_valid;
              }
              const std::size_t MIN_VALID_CELLS = 50;
              if (num_valid < MIN_VALID_CELLS) {
                continue;
              }

              // IRLS profile fit on the grid. `ref` is already normalised
              // over its valid cells during learning and must not be
              // re-normalised here.
              ProfileFitter<double> fit(data_grid.const_ref(),
                                        bg_grid.const_ref(),
                                        combined_mask.const_ref(),
                                        ref,
                                        1e-3,
                                        100);

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

    TransformSpec spec_;
  };

}}  // namespace dials::algorithms

#endif  // DIALS_ALGORITHMS_PROFILE_MODEL_GAUSSIAN_RS_MODELLER_H
