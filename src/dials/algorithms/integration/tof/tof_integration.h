#ifndef DIALS_ALGORITHMS_INTEGRATION_TOF_INTEGRATION_H
#define DIALS_ALGORITHMS_INTEGRATION_TOF_INTEGRATION_H

#include <memory>
#include <boost/optional.hpp>
#include <dxtbx/imageset.h>
#include <dxtbx/format/image.h>
#include <dxtbx/array_family/flex_table.h>
#include <dials/model/data/shoebox.h>
#include <dxtbx/model/detector.h>
#include <dxtbx/model/beam.h>
#include <dxtbx/model/scan.h>
#include <dxtbx/model/goniometer.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/algorithms/integration/processor.h>
#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include <scitbx/constants.h>
#include <dials/model/data/mask_code.h>
#include <dials/algorithms/integration/tof/tof_profile_1d_ibix.h>
#include <dials/algorithms/integration/tof/tof_profile_3d_gutmann.h>
#include <dials/algorithms/scaling/tof/tof_scaling.h>
#include <dials/util/thread_pool.h>

namespace dials { namespace algorithms {

  using dials::algorithms::Shoebox;
  using dials::algorithms::ShoeboxProcessor;
  using dials::model::Background;
  using dials::model::BackgroundUsed;
  using dials::model::Foreground;
  using dials::model::Overlapped;
  using dials::model::Valid;
  using dxtbx::ImageSequence;
  using dxtbx::af::flex_table;
  using dxtbx::model::Detector;
  using dxtbx::model::Experiment;
  using dxtbx::model::Goniometer;
  using dxtbx::model::PolychromaticBeam;
  using dxtbx::model::Scan;
  using dxtbx::model::scan_property_types;
  using scitbx::deg_as_rad;
  using scitbx::mat3;
  using scitbx::vec2;
  using scitbx::vec3;
  using scitbx::af::int6;
  using scitbx::constants::m_n;
  using scitbx::constants::pi;
  using scitbx::constants::Planck;

  // Holds information common to different shoeboxes
  struct TOFGeometryContext {
    Detector detector;
    vec3<double> unit_s0;
    double sample_to_source_distance;
    scitbx::af::shared<double> img_tof;
    vec2<std::size_t> image_size;
    int n_panels;

    TOFGeometryContext(Experiment& experiment, ImageSequence& data)
        : detector(*experiment.get_detector()) {
      Scan scan = *experiment.get_scan();

      std::shared_ptr<dxtbx::model::BeamBase> beam_ptr = experiment.get_beam();
      std::shared_ptr<PolychromaticBeam> beam =
        std::dynamic_pointer_cast<PolychromaticBeam>(beam_ptr);
      DIALS_ASSERT(beam != nullptr);
      unit_s0 = beam->get_unit_s0();
      sample_to_source_distance = beam->get_sample_to_source_distance();

      img_tof = scan.get_property<double>("time_of_flight");

      n_panels = detector.size();
      int num_images = data.size();
      image_size = detector[0].get_image_size();
      DIALS_ASSERT(num_images == img_tof.size());
    }
  };

  // First pass over a shoebox counts foreground/background pixels
  struct ShoeboxPixelCount {
    bool success = true;
    int n_signal = 0;
    int n_background = 0;
    scitbx::af::shared<double> incident_spectrum;
    scitbx::af::shared<double> empty_spectrum;
    scitbx::af::shared<std::size_t> n_contrib;
  };

  inline ShoeboxPixelCount scan_shoebox_pixels(
    const Shoebox<>& shoebox,
    const int6& bbox,
    const vec2<std::size_t>& image_size,
    const Shoebox<>* incident_shoebox = nullptr,
    const Shoebox<>* empty_shoebox = nullptr) {
    ShoeboxPixelCount result;
    bool has_incident = incident_shoebox != nullptr;
    if (has_incident) {
      result.incident_spectrum = scitbx::af::shared<double>(shoebox.zsize(), 0.0);
      result.empty_spectrum = scitbx::af::shared<double>(shoebox.zsize(), 0.0);
      result.n_contrib = scitbx::af::shared<std::size_t>(shoebox.zsize(), 0);
    }

    int bg_code = Valid | Background | BackgroundUsed;

    // Shoebox data are ordered (z, y, x)
    for (std::size_t z = 0; z < shoebox.zsize(); ++z) {
      if (!result.success) {
        break;
      }
      for (std::size_t y = 0; y < shoebox.ysize(); ++y) {
        int panel_y = bbox[2] + y;
        if (panel_y > image_size[1] || panel_y < 0) {
          continue;
        }
        for (std::size_t x = 0; x < shoebox.xsize(); ++x) {
          int panel_x = bbox[0] + x;
          if (panel_x > image_size[0] || panel_x < 0) {
            continue;
          }

          if (has_incident) {
            result.incident_spectrum[z] += incident_shoebox->data(z, y, x);
            result.empty_spectrum[z] += empty_shoebox->data(z, y, x);
            result.n_contrib[z]++;
          }

          int mask = shoebox.mask(z, y, x);
          if ((mask & Foreground) == Foreground) {
            if ((mask & Valid) == Valid && (mask & Overlapped) == 0) {
              result.n_signal++;
            } else {
              result.success = false;
            }
          } else if ((mask & bg_code) == bg_code) {
            result.n_background++;
          }
        }
      }
    }
    return result;
  }

  // Loads imagesequence data into reflection shoeboxes
  inline dials::af::reflection_table load_shoebox_image_data(
    dials::af::reflection_table& reflection_table,
    std::shared_ptr<ImageSequence> reference_data,
    int n_panels,
    int num_images) {
    boost::python::dict d;
    dials::af::reflection_table ref_table =
      dxtbx::af::flex_table_suite::deepcopy(reflection_table, d);

    dials::af::ref<Shoebox<>> ref_shoeboxes = ref_table["shoebox"];
    for (std::size_t i = 0; i < ref_table.size(); ++i) {
      ref_shoeboxes[i].deallocate();
    }

    ShoeboxProcessor processor(ref_table, n_panels, 0, num_images, false);

    for (std::size_t img_num = 0; img_num < static_cast<std::size_t>(num_images);
         ++img_num) {
      dxtbx::format::Image<double> img = reference_data->get_corrected_data(img_num);
      dxtbx::format::Image<bool> mask = reference_data->get_mask(img_num);

      dials::af::shared<scitbx::af::versa<double, scitbx::af::c_grid<2>>> output_data(
        n_panels);
      dials::af::shared<scitbx::af::versa<bool, scitbx::af::c_grid<2>>> output_mask(
        n_panels);

      for (std::size_t i = 0; i < output_data.size(); ++i) {
        output_data[i] = img.tile(i).data();
        output_mask[i] = mask.tile(i).data();
      }
      processor.next_data_only(
        dials::model::Image<double>(output_data.const_ref(), output_mask.const_ref()));
    }

    return ref_table;
  }

  struct PixelCorrection {
    double I = 0;
    double B = 0;
    double var_I = 0;
    double var_B = 0;
    double I_raw = 0;
  };

  class PixelCorrector {
  public:
    PixelCorrector(const TOFGeometryContext& geometry,
                   bool apply_lorentz_correction,
                   int n_signal,
                   int n_background,
                   const dials_scaling::TOFIncidentSpectrumParams* incident_params,
                   const dials_scaling::TOFAbsorptionParams* absorption_params,
                   const scitbx::af::shared<double>* smoothed_incident,
                   const scitbx::af::shared<double>* smoothed_empty,
                   const scitbx::af::shared<std::size_t>* n_contrib)
        : geometry_(geometry),
          apply_lorentz_correction_(apply_lorentz_correction),
          n_signal_(n_signal),
          n_background_(n_background),
          incident_params_(incident_params),
          absorption_params_(absorption_params),
          smoothed_incident_(smoothed_incident),
          smoothed_empty_(smoothed_empty),
          n_contrib_(n_contrib) {}

    // Returns false if the pixel cannot be corrected and should be skipped
    bool correct(int panel,
                 int panel_x,
                 int panel_y,
                 int frame_z,
                 std::size_t z,
                 double tof_us,
                 double raw_S,
                 double raw_B,
                 PixelCorrection& out) const {
      double tof_s = tof_us * std::pow(10, -6);

      double var_S = std::abs(raw_S);
      double var_B = std::abs(raw_B);
      if (n_background_ > 0) {
        var_B *= (1.0 + double(n_signal_) / double(n_background_));
      }

      bool has_incident = incident_params_ != nullptr;
      bool has_absorption = absorption_params_ != nullptr;

      // wl/two_theta are needed for the Lorentz correction and for the
      // absorption correction, so compute them once if either is active.
      double wl = 0.0, two_theta = 0.0, L = 1.0;
      if (apply_lorentz_correction_ || has_absorption) {
        scitbx::vec3<double> s1 = geometry_.detector[panel].get_pixel_lab_coord(
          scitbx::vec2<double>(panel_x, panel_y));
        double distance =
          (s1.length() + geometry_.sample_to_source_distance) * std::pow(10, -3);
        wl = ((Planck * tof_s) / (m_n * distance)) * std::pow(10, 10);
        two_theta = geometry_.detector[panel].get_two_theta_at_pixel(
          geometry_.unit_s0, scitbx::vec2<double>(panel_x, panel_y));
        if (apply_lorentz_correction_) {
          double sin_two_theta_sq = std::pow(sin(two_theta * .5), 2);
          L = sin_two_theta_sq / std::pow(wl, 4);
        }
      }

      // Spherical absorption correction
      double T = 1.0;
      int two_theta_idx = 0;
      if (has_absorption) {
        double two_theta_deg = two_theta * (180.0 / pi);
        two_theta_idx = static_cast<int>(two_theta_deg / 10);
        double sample_muR =
          (absorption_params_->sample_linear_scattering_c
           + (absorption_params_->sample_linear_absorption_c / 1.8) * wl)
          * absorption_params_->sample_radius;
        T = dials_scaling::tof_pixel_spherical_absorption_correction(
          sample_muR, two_theta, two_theta_idx);
        if (T < 1e-7) {
          return false;
        }
      }

      double I0, var_I0, I_raw, B, var_B_pre;
      double J = 1.0, var_J = 0.0;

      // Bin-width normalisation is used when there is no incident/empty
      if (!has_incident) {
        double bin_width = geometry_.img_tof[frame_z] - geometry_.img_tof[frame_z - 1];
        I0 = (raw_S - raw_B) / bin_width;
        var_I0 = (var_S + var_B) / (bin_width * bin_width);
        I_raw = raw_S / bin_width;
        B = raw_B / bin_width;
        var_B_pre = var_B / (bin_width * bin_width);
      } else {
        double n = double((*n_contrib_)[z]);
        double raw_J = (*smoothed_incident_)[z] / n;
        double raw_E = (*smoothed_empty_)[z] / n;

        double sample_pc = incident_params_->sample_proton_charge;
        double incident_pc = incident_params_->incident_proton_charge;
        double empty_pc = incident_params_->empty_proton_charge;

        double S = raw_S / sample_pc;
        double Bn = raw_B / sample_pc;
        double Jn = raw_J / incident_pc;
        double En = raw_E / empty_pc;

        double var_Sn = var_S / (sample_pc * sample_pc);
        double var_Bn = var_B / (sample_pc * sample_pc);
        double var_Jn = (std::abs(raw_J) * n) / (incident_pc * incident_pc * n * n);
        double var_En = (std::abs(raw_E) * n) / (empty_pc * empty_pc * n * n);

        if (has_absorption) {
          // Spherical absorption correction for the incident spectrum
          double incident_muR =
            (absorption_params_->incident_linear_scattering_c
             + (absorption_params_->incident_linear_absorption_c / 1.8) * wl)
            * absorption_params_->incident_radius;
          double J_T = dials_scaling::tof_pixel_spherical_absorption_correction(
            incident_muR, two_theta, two_theta_idx);
          if (J_T < 1e-7) {
            return false;
          }
          Jn /= J_T;
          if (Jn < 1e-7) {
            return false;
          }
          var_Jn /= (J_T * J_T);
        }

        I0 = S - Bn - En;
        var_I0 = var_Sn + var_Bn + var_En;
        I_raw = S;
        B = Bn;
        var_B_pre = var_Bn;
        J = Jn;
        var_J = var_Jn;
      }

      // Apply the Lorentz correction and normalise
      double JT = J * T;
      out.I = L * I0 / JT;
      out.B = L * B / JT;
      out.I_raw = L * I_raw / JT;
      out.var_I =
        (L * L / (JT * JT)) * var_I0 + (L * L * I0 * I0 / std::pow(JT, 4)) * var_J;
      out.var_B =
        (L * L / (JT * JT)) * var_B_pre + (L * L * B * B / std::pow(JT, 4)) * var_J;
      return true;
    }

  private:
    const TOFGeometryContext& geometry_;
    bool apply_lorentz_correction_;
    int n_signal_;
    int n_background_;
    const dials_scaling::TOFIncidentSpectrumParams* incident_params_;
    const dials_scaling::TOFAbsorptionParams* absorption_params_;
    const scitbx::af::shared<double>* smoothed_incident_;
    const scitbx::af::shared<double>* smoothed_empty_;
    const scitbx::af::shared<std::size_t>* n_contrib_;
  };

  struct ShoeboxIntegrationResult {
    bool success = true;
    double intensity = 0.0;
    double variance = 0.0;
    scitbx::af::shared<double> tof_z;
    scitbx::af::shared<double> projected_intensity;
    scitbx::af::shared<double> raw_projected_intensity;
    scitbx::af::shared<double> projected_background;
    scitbx::af::versa<double, af::c_grid<3>> intensity_3d;
    scitbx::af::versa<double, af::c_grid<3>> background_var_3d;
    scitbx::af::versa<vec3<double>, af::c_grid<3>> coords_3d;
  };

  // Second pass over a shoebox
  // applies the pixel corrector and accumulates summation intensity/variance
  // plus every array needed by 1D/3D profile fitting
  inline ShoeboxIntegrationResult integrate_shoebox(const Shoebox<>& shoebox,
                                                    const int6& bbox,
                                                    const TOFGeometryContext& geometry,
                                                    const PixelCorrector& corrector,
                                                    bool initial_success = true) {
    ShoeboxIntegrationResult result;
    result.success = initial_success;
    result.tof_z = scitbx::af::shared<double>(shoebox.zsize(), 0.0);
    result.projected_intensity = scitbx::af::shared<double>(shoebox.zsize(), 0.0);
    result.raw_projected_intensity = scitbx::af::shared<double>(shoebox.zsize(), 0.0);
    result.projected_background = scitbx::af::shared<double>(shoebox.zsize(), 0.0);

    auto acc = shoebox.data.accessor();
    af::c_grid<3> transposed_acc(acc[2], acc[1], acc[0]);
    result.intensity_3d = scitbx::af::versa<double, af::c_grid<3>>(transposed_acc);
    result.background_var_3d = scitbx::af::versa<double, af::c_grid<3>>(transposed_acc);
    result.coords_3d = scitbx::af::versa<vec3<double>, af::c_grid<3>>(transposed_acc);

    int panel = shoebox.panel;

    // Shoebox data are ordered (z, y, x)
    for (std::size_t z = 0; z < shoebox.zsize(); ++z) {
      double intensity_z = 0;
      double intensity_raw_z = 0;
      double background_z = 0;

      if (!result.success) {
        break;
      }

      int frame_z = bbox[4] + z;
      double tof = geometry.img_tof[frame_z];
      result.tof_z[z] = tof;  // (us)

      for (std::size_t y = 0; y < shoebox.ysize(); ++y) {
        int panel_y = bbox[2] + y;
        if (panel_y > geometry.image_size[1] || panel_y < 0) {
          continue;
        }
        for (std::size_t x = 0; x < shoebox.xsize(); ++x) {
          int panel_x = bbox[0] + x;
          if (panel_x > geometry.image_size[0] || panel_x < 0) {
            continue;
          }

          double raw_S = shoebox.data(z, y, x);
          double raw_B = shoebox.background(z, y, x);
          int mask = shoebox.mask(z, y, x);

          PixelCorrection pixel;
          if (!corrector.correct(
                panel, panel_x, panel_y, frame_z, z, tof, raw_S, raw_B, pixel)) {
            continue;
          }

          intensity_z += pixel.I;
          intensity_raw_z += pixel.I_raw;
          background_z += pixel.B;

          result.intensity_3d(x, y, z) = pixel.I;
          result.background_var_3d(x, y, z) = pixel.var_B;
          double x_c = x + shoebox.xoffset() + 0.5;
          double y_c = y + shoebox.yoffset() + 0.5;
          double z_c = z + shoebox.zoffset() + 0.5;
          result.coords_3d(x, y, z) = vec3<double>(x_c, y_c, z_c);

          // Accumulate summation values if pixel in foreground and valid
          if ((mask & Foreground) == Foreground && (mask & Valid) == Valid
              && (mask & Overlapped) == 0) {
            result.intensity += pixel.I;
            result.variance += pixel.var_I;
          } else if ((mask & Foreground) == Foreground) {
            result.success = false;
            break;
          }
        }
      }

      result.projected_intensity[z] = intensity_z;
      result.raw_projected_intensity[z] = intensity_raw_z;
      result.projected_background[z] = background_z;
    }

    return result;
  }

  // Profile fitting strategy
  class ProfileFitter {
  public:
    virtual ~ProfileFitter() {}
    virtual bool fit(const ShoeboxIntegrationResult& shoebox_result,
                     double& I_prf,
                     double& var_prf) = 0;
  };

  class Profile1DIBIXFitter : public ProfileFitter {
  public:
    explicit Profile1DIBIXFitter(TOFProfile1DIBIXParams params) : params_(params) {}

    bool fit(const ShoeboxIntegrationResult& shoebox_result,
             double& I_prf,
             double& var_prf) override {
      var_prf = shoebox_result.variance;  // Use summation variance as approximation
      return fit_profile_1d_ibix(shoebox_result.projected_intensity.const_ref(),
                                 shoebox_result.tof_z.const_ref(),
                                 params_,
                                 I_prf);
    }

  private:
    TOFProfile1DIBIXParams params_;
  };

  class Profile3DGutmannFitter : public ProfileFitter {
  public:
    explicit Profile3DGutmannFitter(TOFProfile3DGutmannParams params)
        : params_(params) {}

    bool fit(const ShoeboxIntegrationResult& shoebox_result,
             double& I_prf,
             double& var_prf) override {
      var_prf = shoebox_result.variance;  // Use summation variance as approximation
      return fit_profile_3d_gutmann(shoebox_result.coords_3d.const_ref(),
                                    shoebox_result.intensity_3d,
                                    shoebox_result.background_var_3d,
                                    params_,
                                    I_prf);
    }

  private:
    TOFProfile3DGutmannParams params_;
  };

  // Runs incident/empty scanning/smoothing for one shoebox and builds
  // the PixelCorrector
  struct ShoeboxCorrectorInputs {
    ShoeboxPixelCount shoebox_pixel_count;
    scitbx::af::shared<double> smoothed_incident;
    scitbx::af::shared<double> smoothed_empty;
  };

  inline ShoeboxCorrectorInputs prepare_shoebox_corrector_inputs(
    const Shoebox<>& shoebox,
    const int6& bbox,
    const vec2<std::size_t>& image_size,
    const Shoebox<>* incident_shoebox,
    const Shoebox<>* empty_shoebox) {
    ShoeboxCorrectorInputs inputs;
    inputs.shoebox_pixel_count =
      scan_shoebox_pixels(shoebox, bbox, image_size, incident_shoebox, empty_shoebox);
    if (incident_shoebox != nullptr) {
      // Smooth incident and empty to avoid dividing by a noisy signal
      inputs.smoothed_empty =
        dials_scaling::savitzky_golay(inputs.shoebox_pixel_count.empty_spectrum, 7, 2);
      inputs.smoothed_incident = dials_scaling::savitzky_golay(
        inputs.shoebox_pixel_count.incident_spectrum, 7, 2);
    }
    return inputs;
  }

  // Updates reflection_table with intensities and variances corrected by
  // any combination of an incident/empty spectrum normalisation, a
  // spherical absorption correction, and a Lorentz correction.
  // profile-fitted intensity is also produced if profile_fitter is
  // given
  inline void integrate_reflection_table(
    dials::af::reflection_table& reflection_table,
    Experiment& experiment,
    ImageSequence& data,
    boost::optional<dials_scaling::TOFIncidentSpectrumParams> incident_params,
    boost::optional<dials_scaling::TOFAbsorptionParams> absorption_params,
    const bool& apply_lorentz_correction,
    int n_threads,
    std::shared_ptr<ProfileFitter> profile_fitter = nullptr) {
    std::size_t n_reflections = reflection_table.size();
    TOFGeometryContext geometry(experiment, data);

    // incident and empty
    dials::af::reflection_table i_reflection_table;
    dials::af::reflection_table e_reflection_table;
    dials::af::shared<Shoebox<>> i_shoeboxes;
    dials::af::shared<Shoebox<>> e_shoeboxes;

    if (incident_params) {
      i_reflection_table = load_shoebox_image_data(reflection_table,
                                                   incident_params->incident_data,
                                                   geometry.n_panels,
                                                   data.size());
      e_reflection_table = load_shoebox_image_data(
        reflection_table, incident_params->empty_data, geometry.n_panels, data.size());
      i_shoeboxes = i_reflection_table["shoebox"];
      e_shoeboxes = e_reflection_table["shoebox"];
    }

    dials::af::shared<Shoebox<>> shoeboxes = reflection_table["shoebox"];
    dials::af::shared<std::size_t> refl_flags =
      reflection_table.get<std::size_t>("flags");
    dials::af::const_ref<int6> bboxes = reflection_table["bbox"];

    // Arrays to store output of integration
    dials::af::shared<bool> succeeded(reflection_table.size());
    dials::af::shared<double> intensities(reflection_table.size());
    dials::af::shared<double> variances(reflection_table.size());
    dials::af::shared<bool> succeeded_prf;
    dials::af::shared<double> intensities_prf;
    dials::af::shared<double> variances_prf;

    if (profile_fitter) {
      succeeded_prf.resize(reflection_table.size());
      intensities_prf.resize(reflection_table.size());
      variances_prf.resize(reflection_table.size());
    }

    auto worker = [&](std::size_t start, std::size_t end) {
      for (std::size_t i = start; i < end; ++i) {
        if (refl_flags[i] & dials::af::DontIntegrate) {
          continue;
        }
        Shoebox<> shoebox = shoeboxes[i];
        int6 bbox = bboxes[i];

        // incident and empty
        Shoebox<> i_shoebox, e_shoebox;
        const Shoebox<>* i_shoebox_ptr = nullptr;
        const Shoebox<>* e_shoebox_ptr = nullptr;
        if (incident_params) {
          i_shoebox = i_shoeboxes[i];
          e_shoebox = e_shoeboxes[i];
          i_shoebox_ptr = &i_shoebox;
          e_shoebox_ptr = &e_shoebox;
        }

        ShoeboxCorrectorInputs inputs = prepare_shoebox_corrector_inputs(
          shoebox, bbox, geometry.image_size, i_shoebox_ptr, e_shoebox_ptr);

        PixelCorrector corrector(
          geometry,
          apply_lorentz_correction,
          inputs.shoebox_pixel_count.n_signal,
          inputs.shoebox_pixel_count.n_background,
          incident_params.get_ptr(),
          absorption_params.get_ptr(),
          incident_params ? &inputs.smoothed_incident : nullptr,
          incident_params ? &inputs.smoothed_empty : nullptr,
          incident_params ? &inputs.shoebox_pixel_count.n_contrib : nullptr);

        ShoeboxIntegrationResult shoebox_result = integrate_shoebox(
          shoebox, bbox, geometry, corrector, inputs.shoebox_pixel_count.success);

        // Overall values for shoebox summation
        succeeded[i] = shoebox_result.success;
        intensities[i] = shoebox_result.intensity;
        variances[i] = shoebox_result.variance;

        if (profile_fitter) {
          bool profile_success = false;
          double I_prf = 0.0, var_prf = 0.0;
          if (shoebox_result.success) {
            profile_success = profile_fitter->fit(shoebox_result, I_prf, var_prf);
          }
          if (profile_success) {
            intensities_prf[i] = I_prf;
            variances_prf[i] = var_prf;
          }
          succeeded_prf[i] = profile_success;
        }
      }
    };

    dials::util::ThreadPool pool(n_threads);
    std::size_t chunk_size = (n_reflections + n_threads - 1) / n_threads;

    for (int t = 0; t < n_threads; ++t) {
      std::size_t start = t * chunk_size;
      std::size_t end = std::min(start + chunk_size, n_reflections);
      if (start >= end) break;

      pool.post([=]() { worker(start, end); });
    }

    pool.wait();

    reflection_table["intensity.sum.value"] = intensities;
    reflection_table["intensity.sum.variance"] = variances;

    if (profile_fitter) {
      reflection_table["intensity.prf.value"] = intensities_prf;
      reflection_table["intensity.prf.variance"] = variances_prf;
    }

    // Update flags
    dials::af::ref<std::size_t> flags =
      reflection_table.get<std::size_t>("flags").ref();
    for (std::size_t i = 0; i < flags.size(); ++i) {
      if (succeeded[i]) {
        flags[i] &= ~dials::af::FailedDuringSummation;
        flags[i] |= dials::af::IntegratedSum;
      } else {
        flags[i] &= ~dials::af::IntegratedSum;
        flags[i] |= dials::af::FailedDuringSummation;
      }
      if (profile_fitter) {
        if (succeeded_prf[i]) {
          flags[i] &= ~dials::af::FailedDuringProfileFitting;
          flags[i] |= dials::af::IntegratedPrf;
        } else {
          flags[i] &= ~dials::af::IntegratedPrf;
          flags[i] |= dials::af::FailedDuringProfileFitting;
        }
      }
    }
  }

  // Single-reflection variant used for diagnostics/plotting
  // Calculates raw_projected_intensity, projected_intensity, projected_background,
  // sum_intensity and sum_variance, with the same correction
  // combinations as integrate_reflection_table
  inline boost::python::tuple calculate_line_profile_for_reflection(
    dials::af::reflection_table& reflection,
    Experiment& experiment,
    ImageSequence& data,
    boost::optional<dials_scaling::TOFIncidentSpectrumParams> incident_params,
    boost::optional<dials_scaling::TOFAbsorptionParams> absorption_params,
    scitbx::af::shared<double> raw_projected_intensity_out,
    scitbx::af::shared<double> projected_intensity_out,
    scitbx::af::shared<double> projected_background_out,
    scitbx::af::shared<double> tof_z_out,
    const bool& apply_lorentz_correction) {
    // This is a slight hack to make using current interfaces easier
    // E.g. the shoebox processor
    DIALS_ASSERT(reflection.size() == 1);

    TOFGeometryContext geometry(experiment, data);

    dials::af::shared<Shoebox<>> shoeboxes = reflection["shoebox"];
    Shoebox<> shoebox = shoeboxes[0];
    int6 bbox = shoebox.bbox;

    DIALS_ASSERT(raw_projected_intensity_out.size() == shoebox.zsize());
    DIALS_ASSERT(projected_intensity_out.size() == shoebox.zsize());
    DIALS_ASSERT(projected_background_out.size() == shoebox.zsize());
    DIALS_ASSERT(tof_z_out.size() == shoebox.zsize());

    dials::af::reflection_table i_reflection_table;
    dials::af::reflection_table e_reflection_table;
    Shoebox<> i_shoebox, e_shoebox;
    const Shoebox<>* i_shoebox_ptr = nullptr;
    const Shoebox<>* e_shoebox_ptr = nullptr;
    if (incident_params) {
      i_reflection_table = load_shoebox_image_data(
        reflection, incident_params->incident_data, geometry.n_panels, data.size());
      e_reflection_table = load_shoebox_image_data(
        reflection, incident_params->empty_data, geometry.n_panels, data.size());
      dials::af::ref<Shoebox<>> i_shoeboxes = i_reflection_table["shoebox"];
      dials::af::ref<Shoebox<>> e_shoeboxes = e_reflection_table["shoebox"];
      i_shoebox = i_shoeboxes[0];
      e_shoebox = e_shoeboxes[0];
      i_shoebox_ptr = &i_shoebox;
      e_shoebox_ptr = &e_shoebox;
    }

    ShoeboxCorrectorInputs inputs = prepare_shoebox_corrector_inputs(
      shoebox, bbox, geometry.image_size, i_shoebox_ptr, e_shoebox_ptr);

    PixelCorrector corrector(
      geometry,
      apply_lorentz_correction,
      inputs.shoebox_pixel_count.n_signal,
      inputs.shoebox_pixel_count.n_background,
      incident_params.get_ptr(),
      absorption_params.get_ptr(),
      incident_params ? &inputs.smoothed_incident : nullptr,
      incident_params ? &inputs.smoothed_empty : nullptr,
      incident_params ? &inputs.shoebox_pixel_count.n_contrib : nullptr);

    ShoeboxIntegrationResult result = integrate_shoebox(
      shoebox, bbox, geometry, corrector, inputs.shoebox_pixel_count.success);

    for (std::size_t z = 0; z < shoebox.zsize(); ++z) {
      raw_projected_intensity_out[z] = result.raw_projected_intensity[z];
      projected_intensity_out[z] = result.projected_intensity[z];
      projected_background_out[z] = result.projected_background[z];
      tof_z_out[z] = result.tof_z[z];
    }

    return boost::python::make_tuple(result.intensity, result.variance, result.success);
  }

  // As above, but also fits a 1D profile and returns the profile-
  // fitted intensity, the line profile itself (line_profile_out), and the
  // summation results.
  inline boost::python::tuple calculate_line_profile_for_reflection(
    dials::af::reflection_table& reflection,
    Experiment& experiment,
    ImageSequence& data,
    boost::optional<dials_scaling::TOFIncidentSpectrumParams> incident_params,
    boost::optional<dials_scaling::TOFAbsorptionParams> absorption_params,
    scitbx::af::shared<double> raw_projected_intensity_out,
    scitbx::af::shared<double> projected_intensity_out,
    scitbx::af::shared<double> projected_background_out,
    scitbx::af::shared<double> tof_z_out,
    scitbx::af::shared<double> line_profile_out,
    const bool& apply_lorentz_correction,
    TOFProfile1DIBIXParams& profile_params_1d_ibix) {
    boost::python::tuple result =
      calculate_line_profile_for_reflection(reflection,
                                            experiment,
                                            data,
                                            incident_params,
                                            absorption_params,
                                            raw_projected_intensity_out,
                                            projected_intensity_out,
                                            projected_background_out,
                                            tof_z_out,
                                            apply_lorentz_correction);

    double I_sum = boost::python::extract<double>(result[0]);
    double var_sum = boost::python::extract<double>(result[1]);
    bool success = boost::python::extract<bool>(result[2]);

    double I_prf = 0;
    double var_prf = 0;
    bool profile_success = false;

    if (success) {
      profile_success = fit_profile_1d_ibix(projected_intensity_out.const_ref(),
                                            tof_z_out.const_ref(),
                                            profile_params_1d_ibix,
                                            I_prf,
                                            line_profile_out,
                                            true);
    }
    return boost::python::make_tuple(I_prf, var_prf, I_sum, var_sum, profile_success);
  }

  // As calculate_line_profile_for_reflection, but also fits a 3D
  // profile and returns the profile-fitted intensity, the fitted 3D
  // profile, and the summation results.
  inline boost::python::tuple calculate_line_profile_for_reflection_3d(
    dials::af::reflection_table& reflection,
    Experiment& experiment,
    ImageSequence& data,
    boost::optional<dials_scaling::TOFIncidentSpectrumParams> incident_params,
    boost::optional<dials_scaling::TOFAbsorptionParams> absorption_params,
    scitbx::af::shared<double> raw_projected_intensity_out,
    scitbx::af::shared<double> projected_intensity_out,
    scitbx::af::shared<double> projected_background_out,
    scitbx::af::shared<double> tof_z_out,
    const bool& apply_lorentz_correction,
    TOFProfile3DGutmannParams& profile_params_3d_gutmann) {
    DIALS_ASSERT(reflection.size() == 1);

    TOFGeometryContext geometry(experiment, data);

    dials::af::shared<Shoebox<>> shoeboxes = reflection["shoebox"];
    Shoebox<> shoebox = shoeboxes[0];
    int6 bbox = shoebox.bbox;

    DIALS_ASSERT(raw_projected_intensity_out.size() == shoebox.zsize());
    DIALS_ASSERT(projected_intensity_out.size() == shoebox.zsize());
    DIALS_ASSERT(projected_background_out.size() == shoebox.zsize());
    DIALS_ASSERT(tof_z_out.size() == shoebox.zsize());

    dials::af::reflection_table i_reflection_table;
    dials::af::reflection_table e_reflection_table;
    Shoebox<> i_shoebox, e_shoebox;
    const Shoebox<>* i_shoebox_ptr = nullptr;
    const Shoebox<>* e_shoebox_ptr = nullptr;
    if (incident_params) {
      i_reflection_table = load_shoebox_image_data(
        reflection, incident_params->incident_data, geometry.n_panels, data.size());
      e_reflection_table = load_shoebox_image_data(
        reflection, incident_params->empty_data, geometry.n_panels, data.size());
      dials::af::ref<Shoebox<>> i_shoeboxes = i_reflection_table["shoebox"];
      dials::af::ref<Shoebox<>> e_shoeboxes = e_reflection_table["shoebox"];
      i_shoebox = i_shoeboxes[0];
      e_shoebox = e_shoeboxes[0];
      i_shoebox_ptr = &i_shoebox;
      e_shoebox_ptr = &e_shoebox;
    }

    ShoeboxCorrectorInputs inputs = prepare_shoebox_corrector_inputs(
      shoebox, bbox, geometry.image_size, i_shoebox_ptr, e_shoebox_ptr);

    PixelCorrector corrector(
      geometry,
      apply_lorentz_correction,
      inputs.shoebox_pixel_count.n_signal,
      inputs.shoebox_pixel_count.n_background,
      incident_params.get_ptr(),
      absorption_params.get_ptr(),
      incident_params ? &inputs.smoothed_incident : nullptr,
      incident_params ? &inputs.smoothed_empty : nullptr,
      incident_params ? &inputs.shoebox_pixel_count.n_contrib : nullptr);

    ShoeboxIntegrationResult result = integrate_shoebox(
      shoebox, bbox, geometry, corrector, inputs.shoebox_pixel_count.success);

    for (std::size_t z = 0; z < shoebox.zsize(); ++z) {
      raw_projected_intensity_out[z] = result.raw_projected_intensity[z];
      projected_intensity_out[z] = result.projected_intensity[z];
      projected_background_out[z] = result.projected_background[z];
      tof_z_out[z] = result.tof_z[z];
    }

    double I_prf = 0;
    double var_prf = 0;
    bool profile_success = false;
    scitbx::af::versa<double, scitbx::af::c_grid<3>> profile_3d_out(
      result.intensity_3d.accessor());

    if (result.success) {
      profile_success = fit_profile_3d_gutmann(result.coords_3d.const_ref(),
                                               result.intensity_3d,
                                               result.background_var_3d,
                                               profile_params_3d_gutmann,
                                               I_prf,
                                               profile_3d_out,
                                               true);
    }
    return boost::python::make_tuple(I_prf,
                                     var_prf,
                                     result.intensity,
                                     result.variance,
                                     profile_success,
                                     profile_3d_out);
  }

}}  // namespace dials::algorithms

#endif /* DIALS_ALGORITHMS_INTEGRATION_TOF_INTEGRATION_H */
