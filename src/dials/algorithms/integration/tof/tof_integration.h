#ifndef DIALS_ALGORITHMS_INTEGRATION_TOF_INTEGRATION_H
#define DIALS_ALGORITHMS_INTEGRATION_TOF_INTEGRATION_H

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
#include <iostream>
#include <vector>
#include <dials/algorithms/integration/tof/tof_profile1d.h>
#include <dials/algorithms/integration/tof/tof_profile3d.h>
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

  void integrate_reflection_table(
    dials::af::reflection_table &reflection_table,
    Experiment &experiment,
    ImageSequence &data,
    const bool &apply_lorentz_correction,
    int n_threads,
    boost::optional<TOFProfile1DParams> profile_params_1d = boost::none,
    boost::optional<TOFProfile3DParams> profile_params_3d = boost::none) {
    /*
     * Updates reflection_table with intensities and variances with
     * optional Lorentz correction
     */

    // Only one profile fitting method allowed
    DIALS_ASSERT(!(profile_params_1d && profile_params_3d));

    std::size_t n_reflections = reflection_table.size();

    Detector detector = *experiment.get_detector();
    Scan scan = *experiment.get_scan();

    // Required beam params
    std::shared_ptr<dxtbx::model::BeamBase> beam_ptr = experiment.get_beam();
    std::shared_ptr<PolychromaticBeam> beam =
      std::dynamic_pointer_cast<PolychromaticBeam>(beam_ptr);
    DIALS_ASSERT(beam != nullptr);
    vec3<double> unit_s0 = beam->get_unit_s0();
    double sample_to_source_distance = beam->get_sample_to_source_distance();

    // Required scan params
    scitbx::af::shared<double> img_tof = scan.get_property<double>("time_of_flight");

    // Required detector params
    int n_panels = detector.size();
    int num_images = data.size();
    vec2<std::size_t> image_size = detector[0].get_image_size();
    DIALS_ASSERT(num_images == img_tof.size());

    dials::af::shared<Shoebox<>> shoeboxes = reflection_table["shoebox"];
    dials::af::const_ref<int6> bboxes = reflection_table["bbox"];

    // Arrays to store ouput of integration
    dials::af::shared<bool> succeeded(reflection_table.size());
    dials::af::shared<double> intensities(reflection_table.size());
    dials::af::shared<double> variances(reflection_table.size());
    dials::af::shared<bool> succeeded_prf;
    dials::af::shared<double> intensities_prf;
    dials::af::shared<double> variances_prf;

    if (profile_params_1d || profile_params_3d) {
      succeeded_prf.resize(reflection_table.size());
      intensities_prf.resize(reflection_table.size());
      variances_prf.resize(reflection_table.size());
    }

    int bg_code = Valid | Background | BackgroundUsed;

    auto worker = [&](std::size_t start, std::size_t end) {
      for (std::size_t i = start; i < end; ++i) {
        Shoebox<> shoebox = shoeboxes[i];
        int6 bbox = bboxes[i];
        int panel = shoebox.panel;

        bool success = true;
        int n_background = 0;
        int n_signal = 0;
        double intensity = 0.0;
        double variance = 0.0;

        // First pass to get, n_signal, n_background
        // Shoebox data are ordered (z, y, x)
        for (std::size_t z = 0; z < shoebox.zsize(); ++z) {
          if (!success) {
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

              int mask = shoebox.mask(z, y, x);
              if ((mask & Foreground) == Foreground) {
                if ((mask & Valid) == Valid && (mask & Overlapped) == 0) {
                  n_signal++;
                } else {
                  success = false;
                }
              } else if ((mask & bg_code) == bg_code) {
                n_background++;
              }
            }
          }
        }

        // 1D profile fitting params
        scitbx::af::shared<double> tof_z(shoebox.zsize(), 0.0);
        scitbx::af::shared<double> projected_intensity(shoebox.zsize(), 0.0);

        // 3D profile fitting params done in xyz
        auto acc = shoebox.data.accessor();
        af::c_grid<3> transposed_acc(acc[2], acc[1], acc[0]);
        scitbx::af::versa<double, af::c_grid<3>> intensity_3d(transposed_acc);
        scitbx::af::versa<double, af::c_grid<3>> background_var_3d(transposed_acc);
        scitbx::af::versa<vec3<double>, af::c_grid<3>> coords_3d(transposed_acc);

        // Second pass to perform actual integration
        for (std::size_t z = 0; z < shoebox.zsize(); ++z) {
          double intensity_z = 0;

          if (!success) {
            break;
          }

          int frame_z = bbox[4] + z;
          double tof = img_tof[frame_z];
          tof_z[z] = tof;           // (μs)
          tof *= std::pow(10, -6);  // (s)

          for (std::size_t y = 0; y < shoebox.ysize(); ++y) {
            int panel_y = bbox[2] + y;

            // Bounds check
            if (panel_y > image_size[1] || panel_y < 0) {
              continue;
            }

            for (std::size_t x = 0; x < shoebox.xsize(); ++x) {
              int panel_x = bbox[0] + x;

              // Bounds check
              if (panel_x > image_size[0] || panel_x < 0) {
                continue;
              }

              // Raw counts
              double raw_S = shoebox.data(z, y, x);
              double raw_B = shoebox.background(z, y, x);
              int mask = shoebox.mask(z, y, x);

              // Variances
              double var_S = std::abs(raw_S);
              double var_B = std::abs(var_B);
              if (n_background > 0) {
                var_B *= (1.0 + double(n_signal) / double(n_background));
              }

              double L;

              // Lorentz correction
              if (apply_lorentz_correction) {
                scitbx::vec3<double> s1 = detector[panel].get_pixel_lab_coord(
                  scitbx::vec2<double>(panel_x, panel_y));
                double distance = s1.length() + sample_to_source_distance;
                distance *= std::pow(10, -3);  // (m)
                double wl = ((Planck * tof) / (m_n * (distance))) * std::pow(10, 10);
                double two_theta = detector[panel].get_two_theta_at_pixel(
                  unit_s0, scitbx::vec2<double>(panel_x, panel_y));
                double sin_two_theta_sq = std::pow(sin(two_theta * .5), 2);
                L = sin_two_theta_sq / std::pow(wl, 4);
              }

              // Account for some data having different sized bin widths
              double bin_width_correction = img_tof[frame_z] - img_tof[frame_z - 1];

              // Net signal
              double I0 = raw_S - raw_B;
              double var_I0 = var_S + var_B;

              double I = I0 / bin_width_correction;
              double B = raw_B / bin_width_correction;
              double var_I = var_I0 / (bin_width_correction * bin_width_correction);
              var_B /= (bin_width_correction * bin_width_correction);

              if (apply_lorentz_correction) {
                I *= L;
                B *= L;
                var_I *= (L * L);
                var_B *= (L * L);
              }

              intensity_z += I;

              if (profile_params_3d) {
                intensity_3d(x, y, z) = I;
                background_var_3d(x, y, z) = var_B;
                double x_c = x + shoebox.xoffset() + 0.5;
                double y_c = y + shoebox.yoffset() + 0.5;
                double z_c = z + shoebox.zoffset() + 0.5;
                coords_3d(x, y, z) = vec3<double>(x_c, y_c, z_c);
              }

              // Accumulate summation values if pixel in foreground and valid
              if ((mask & Foreground) == Foreground && (mask & Valid) == Valid
                  && (mask & Overlapped) == 0) {
                intensity += I;
                variance += var_I;
              } else if ((mask & Foreground) == Foreground) {
                success = false;
                break;
              }
            }
          }

          if (profile_params_1d) {  // Doing 1D profile fitting
            projected_intensity[z] = intensity_z;
          }
        }

        // Overall values for shoebox summation
        succeeded[i] = success;
        intensities[i] = intensity;
        variances[i] = variance;

        if (profile_params_1d && success) {
          bool profile_success = false;
          double I_prf;
          profile_success = fit_profile1d(projected_intensity.const_ref(),
                                          tof_z.const_ref(),
                                          *profile_params_1d,
                                          I_prf);
          if (profile_success) {
            intensities_prf[i] = I_prf;
            variances_prf[i] = variance;  // Use summation variance as approximation
          }
          succeeded_prf[i] = profile_success;

        } else if (profile_params_3d && success) {
          bool profile_success = false;

          double I_prf, var_prf;
          profile_success = fit_profile3d(coords_3d.const_ref(),
                                          intensity_3d,
                                          background_var_3d,
                                          *profile_params_3d,
                                          I_prf);
          if (profile_success) {
            intensities_prf[i] = I_prf;
            variances_prf[i] = variance;  // Use summation variance as approximation
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

    if (profile_params_1d || profile_params_3d) {
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
      if (profile_params_1d || profile_params_3d) {
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

  void integrate_reflection_table(
    dials::af::reflection_table &reflection_table,
    Experiment &experiment,
    ImageSequence &data,
    const dials_scaling::TOFIncidentSpectrumParams &incident_params,
    const bool &apply_lorentz_correction,
    int n_threads,
    boost::optional<TOFProfile1DParams> profile_params_1d = boost::none,
    boost::optional<TOFProfile3DParams> profile_params_3d = boost::none) {
    /*
     * Updates reflection_table with intensities and variances corrected by
     * incident and empty runs, and an optional Lorentz correction
     */

    // Only one profile fitting method allowed
    DIALS_ASSERT(!(profile_params_1d && profile_params_3d));

    std::size_t n_reflections = reflection_table.size();
    Detector detector = *experiment.get_detector();
    Scan scan = *experiment.get_scan();

    // Required beam params
    std::shared_ptr<dxtbx::model::BeamBase> beam_ptr = experiment.get_beam();
    std::shared_ptr<PolychromaticBeam> beam =
      std::dynamic_pointer_cast<PolychromaticBeam>(beam_ptr);
    DIALS_ASSERT(beam != nullptr);
    vec3<double> unit_s0 = beam->get_unit_s0();
    double sample_to_source_distance = beam->get_sample_to_source_distance();

    // Required scan params
    scitbx::af::shared<double> img_tof = scan.get_property<double>("time_of_flight");

    // Required detector params
    int n_panels = detector.size();
    int num_images = data.size();
    vec2<std::size_t> image_size = detector[0].get_image_size();
    DIALS_ASSERT(num_images == img_tof.size());

    // Copy reflections for incident and empty runs
    boost::python::dict d;
    dials::af::reflection_table i_reflection_table =
      dxtbx::af::flex_table_suite::deepcopy(reflection_table, d);
    dials::af::reflection_table e_reflection_table =
      dxtbx::af::flex_table_suite::deepcopy(reflection_table, d);

    dials::af::ref<Shoebox<>> i_shoeboxes = i_reflection_table["shoebox"];
    dials::af::ref<Shoebox<>> e_shoeboxes = e_reflection_table["shoebox"];
    for (std::size_t i = 0; i < i_reflection_table.size(); ++i) {
      i_shoeboxes[i].deallocate();
      e_shoeboxes[i].deallocate();
    }

    ShoeboxProcessor incident_shoebox_processor(
      i_reflection_table, n_panels, 0, num_images, false);

    ShoeboxProcessor empty_shoebox_processor(
      e_reflection_table, n_panels, 0, num_images, false);

    // Get shoeboxes for incident data
    for (std::size_t img_num = 0; img_num < num_images; ++img_num) {
      dxtbx::format::Image<double> img =
        incident_params.incident_data->get_corrected_data(img_num);
      dxtbx::format::Image<bool> mask =
        incident_params.incident_data->get_mask(img_num);

      dials::af::shared<scitbx::af::versa<double, scitbx::af::c_grid<2>>> output_data(
        n_panels);
      dials::af::shared<scitbx::af::versa<bool, scitbx::af::c_grid<2>>> output_mask(
        n_panels);

      for (std::size_t i = 0; i < output_data.size(); ++i) {
        output_data[i] = img.tile(i).data();
        output_mask[i] = mask.tile(i).data();
      }
      incident_shoebox_processor.next_data_only(
        dials::model::Image<double>(output_data.const_ref(), output_mask.const_ref()));
    }

    // Get shoeboxes for empty data
    for (std::size_t img_num = 0; img_num < num_images; ++img_num) {
      dxtbx::format::Image<double> img =
        incident_params.empty_data->get_corrected_data(img_num);
      dxtbx::format::Image<bool> mask = incident_params.empty_data->get_mask(img_num);

      dials::af::shared<scitbx::af::versa<double, scitbx::af::c_grid<2>>> output_data(
        n_panels);
      dials::af::shared<scitbx::af::versa<bool, scitbx::af::c_grid<2>>> output_mask(
        n_panels);

      for (std::size_t i = 0; i < output_data.size(); ++i) {
        output_data[i] = img.tile(i).data();
        output_mask[i] = mask.tile(i).data();
      }
      empty_shoebox_processor.next_data_only(
        dials::model::Image<double>(output_data.const_ref(), output_mask.const_ref()));
    }

    dials::af::shared<Shoebox<>> shoeboxes = reflection_table["shoebox"];
    dials::af::const_ref<int6> bboxes = reflection_table["bbox"];

    // Arrays to store ouput of integration
    dials::af::shared<bool> succeeded(reflection_table.size());
    dials::af::shared<double> intensities(reflection_table.size());
    dials::af::shared<double> variances(reflection_table.size());
    dials::af::shared<bool> succeeded_prf;
    dials::af::shared<double> intensities_prf;
    dials::af::shared<double> variances_prf;

    if (profile_params_1d || profile_params_3d) {
      succeeded_prf.resize(reflection_table.size());
      intensities_prf.resize(reflection_table.size());
      variances_prf.resize(reflection_table.size());
    }

    int bg_code = Valid | Background | BackgroundUsed;

    auto worker = [&](std::size_t start, std::size_t end) {
      for (std::size_t i = start; i < end; ++i) {
        Shoebox<> shoebox = shoeboxes[i];
        Shoebox<> i_shoebox = i_shoeboxes[i];
        Shoebox<> e_shoebox = e_shoeboxes[i];
        int6 bbox = bboxes[i];
        int panel = shoebox.panel;

        bool success = true;
        int n_background = 0;
        int n_signal = 0;
        double intensity = 0.0;
        double variance = 0.0;

        // First pass to get, n_signal, n_background, incident spectrum, empty spectrum
        // Shoebox data are ordered (z, y, x)
        scitbx::af::shared<double> incident_spectrum(shoebox.zsize(), 0.0);
        scitbx::af::shared<double> empty_spectrum(shoebox.zsize(), 0.0);
        std::vector<std::size_t> n_contrib(shoebox.zsize(), 0);

        for (std::size_t z = 0; z < shoebox.zsize(); ++z) {
          if (!success) {
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

              incident_spectrum[z] += i_shoebox.data(z, y, x);
              empty_spectrum[z] += e_shoebox.data(z, y, x);
              n_contrib[z]++;

              int mask = shoebox.mask(z, y, x);
              if ((mask & Foreground) == Foreground) {
                if ((mask & Valid) == Valid && (mask & Overlapped) == 0) {
                  n_signal++;
                } else {
                  success = false;
                }
              } else if ((mask & bg_code) == bg_code) {
                n_background++;
              }
            }
          }
        }

        // Smooth incident and empty to avoid dividing by a noisy signal
        scitbx::af::shared<double> smoothed_empty =
          dials_scaling::savitzky_golay(empty_spectrum, 7, 2);
        scitbx::af::shared<double> smoothed_incident =
          dials_scaling::savitzky_golay(incident_spectrum, 7, 2);

        // 1D profile fitting params
        scitbx::af::shared<double> tof_z(shoebox.zsize(), 0.0);
        scitbx::af::shared<double> projected_intensity(shoebox.zsize(), 0.0);

        // 3D profile fitting params done in xyz
        auto acc = shoebox.data.accessor();
        af::c_grid<3> transposed_acc(acc[2], acc[1], acc[0]);
        scitbx::af::versa<double, af::c_grid<3>> intensity_3d(transposed_acc);
        scitbx::af::versa<double, af::c_grid<3>> background_var_3d(transposed_acc);
        scitbx::af::versa<vec3<double>, af::c_grid<3>> coords_3d(transposed_acc);

        // Second pass to perform actual integration
        for (std::size_t z = 0; z < shoebox.zsize(); ++z) {
          double intensity_z = 0;

          if (!success) {
            break;
          }

          int frame_z = bbox[4] + z;
          double tof = img_tof[frame_z];
          tof_z[z] = tof;
          tof *= std::pow(10, -6);  // (s)

          for (std::size_t y = 0; y < shoebox.ysize(); ++y) {
            int panel_y = bbox[2] + y;
            // Bounds check
            if (panel_y > image_size[1] || panel_y < 0) {
              continue;
            }
            for (std::size_t x = 0; x < shoebox.xsize(); ++x) {
              int panel_x = bbox[0] + x;
              // Bounds check
              if (panel_x > image_size[0] || panel_x < 0) {
                continue;
              }

              // Raw counts
              double raw_S = shoebox.data(z, y, x);
              double raw_B = shoebox.background(z, y, x);
              double n = n_contrib[z];
              double raw_J = smoothed_incident[z] / n;
              double raw_E = smoothed_empty[z] / n;
              int mask = shoebox.mask(z, y, x);

              // Normalize by proton charge
              double sample_pc = incident_params.sample_proton_charge;
              double incident_pc = incident_params.incident_proton_charge;
              double empty_pc = incident_params.empty_proton_charge;
              double S = raw_S / sample_pc;
              double B = raw_B / sample_pc;
              double J = raw_J / incident_pc;
              double E = raw_E / empty_pc;

              // Variances
              double var_S = std::abs(raw_S) / (sample_pc * sample_pc);
              double var_B = std::abs(raw_B) / (sample_pc * sample_pc);
              if (n_background > 0) {
                var_B *= (1.0 + double(n_signal) / double(n_background));
              }
              double var_J =
                (std::abs(raw_J) * n) / (incident_pc * incident_pc * n * n);
              double var_E = (std::abs(raw_E) * n) / (empty_pc * empty_pc * n * n);

              // Lorentz correction
              double L;
              if (apply_lorentz_correction) {
                scitbx::vec3<double> s1 = detector[panel].get_pixel_lab_coord(
                  scitbx::vec2<double>(panel_x, panel_y));
                double distance = s1.length() + sample_to_source_distance;
                distance *= std::pow(10, -3);  // (m)
                double wl = ((Planck * tof) / (m_n * (distance))) * std::pow(10, 10);
                double two_theta = detector[panel].get_two_theta_at_pixel(
                  unit_s0, scitbx::vec2<double>(panel_x, panel_y));
                double sin_two_theta_sq = std::pow(sin(two_theta * .5), 2);
                L = sin_two_theta_sq / std::pow(wl, 4);
              }

              // Net signal
              double I0 = S - B - E;
              double var_I0 = var_S + var_B + var_E;

              // Apply Lorentz correction and normalise by incident spectrum
              // (bin_width correction not needed as this cancels with incident spectrum
              // division)
              double I, var_I;
              if (apply_lorentz_correction) {
                I = L * I0 / J;
                B *= L / J;
                // Var( L*I0/J ) = (L/J)^2 Var(I0) + (L*I0/J^2)^2 Var(J)
                var_I = (L * L / (J * J)) * var_I0
                        + (L * L * I0 * I0 / (J * J * J * J)) * var_J;
                var_B =
                  (L * L / (J * J)) * var_B + (L * L * B * B / (J * J * J * J)) * var_J;
              } else {
                I = I0 / J;
                B /= J;
                // Var( I0/J ) = (1/J)^2 Var(I0) + (I0/J^2)^2 Var(J)
                var_I = (1 / (J * J)) * var_I0 + (I0 * I0 / (J * J * J * J)) * var_J;
                var_B = (1 / (J * J)) * var_B + (B * B / (J * J * J * J)) * var_J;
              }

              intensity_z += I;

              if (profile_params_3d) {
                intensity_3d(x, y, z) = I;
                background_var_3d(x, y, z) = var_B;
                double x_c = x + shoebox.xoffset() + 0.5;
                double y_c = y + shoebox.yoffset() + 0.5;
                double z_c = z + shoebox.zoffset() + 0.5;
                coords_3d(x, y, z) = vec3<double>(x_c, y_c, z_c);
              }

              // Accumulate summation values if pixel in foreground and valid
              if ((mask & Foreground) == Foreground && (mask & Valid) == Valid
                  && (mask & Overlapped) == 0) {
                intensity += I;
                variance += var_I;
              } else if ((mask & Foreground) == Foreground) {
                success = false;
                break;
              }
            }
          }
          if (profile_params_1d) {
            projected_intensity[z] = intensity_z;
          }
        }

        // Overall values for shoebox summation
        succeeded[i] = success;
        intensities[i] = intensity;
        variances[i] = variance;

        if (profile_params_1d) {
          bool profile_success = false;
          double I_prf;
          profile_success = fit_profile1d(projected_intensity.const_ref(),
                                          tof_z.const_ref(),
                                          *profile_params_1d,
                                          I_prf);
          if (profile_success) {
            intensities_prf[i] = I_prf;
            variances_prf[i] = variance;  // Use summation variance as approximation
          }
          succeeded_prf[i] = profile_success;
        } else if (profile_params_3d) {
          bool profile_success = false;

          double I_prf, var_prf;
          profile_success = fit_profile3d(coords_3d.const_ref(),
                                          intensity_3d,
                                          background_var_3d,
                                          *profile_params_3d,
                                          I_prf);
          if (profile_success) {
            intensities_prf[i] = I_prf;
            variances_prf[i] = variance;  // Use summation variance as approximation
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

    if (profile_params_1d || profile_params_3d) {
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
      if (profile_params_1d || profile_params_3d) {
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

  void integrate_reflection_table(
    dials::af::reflection_table &reflection_table,
    Experiment &experiment,
    ImageSequence &data,
    const dials_scaling::TOFIncidentSpectrumParams &incident_params,
    const dials_scaling::TOFAbsorptionParams &corrections_data,
    const bool &apply_lorentz_correction,
    int n_threads,
    boost::optional<TOFProfile1DParams> profile_params_1d = boost::none,
    boost::optional<TOFProfile3DParams> profile_params_3d = boost::none) {
    /*
     * Updates reflection_table with intensities and variances corrected by
     * incident and empty runs, a spherical absorption correction,
     * and an optional Lorentz correction
     */

    // Only one profile fitting method allowed
    DIALS_ASSERT(!(profile_params_1d && profile_params_3d));

    std::size_t n_reflections = reflection_table.size();

    Detector detector = *experiment.get_detector();
    Scan scan = *experiment.get_scan();

    // Required beam params
    std::shared_ptr<dxtbx::model::BeamBase> beam_ptr = experiment.get_beam();
    std::shared_ptr<PolychromaticBeam> beam =
      std::dynamic_pointer_cast<PolychromaticBeam>(beam_ptr);
    DIALS_ASSERT(beam != nullptr);
    vec3<double> unit_s0 = beam->get_unit_s0();
    double sample_to_source_distance = beam->get_sample_to_source_distance();

    // Required scan params
    scitbx::af::shared<double> img_tof = scan.get_property<double>("time_of_flight");

    // Required detector params
    int n_panels = detector.size();
    int num_images = data.size();
    vec2<std::size_t> image_size = detector[0].get_image_size();
    DIALS_ASSERT(num_images == img_tof.size());

    // Copy reflections for incident and empty runs
    boost::python::dict d;
    dials::af::reflection_table i_reflection_table =
      dxtbx::af::flex_table_suite::deepcopy(reflection_table, d);
    dials::af::reflection_table e_reflection_table =
      dxtbx::af::flex_table_suite::deepcopy(reflection_table, d);

    dials::af::ref<Shoebox<>> i_shoeboxes = i_reflection_table["shoebox"];
    dials::af::ref<Shoebox<>> e_shoeboxes = e_reflection_table["shoebox"];
    for (std::size_t i = 0; i < i_reflection_table.size(); ++i) {
      i_shoeboxes[i].deallocate();
      e_shoeboxes[i].deallocate();
    }

    ShoeboxProcessor incident_shoebox_processor(
      i_reflection_table, n_panels, 0, num_images, false);

    ShoeboxProcessor empty_shoebox_processor(
      e_reflection_table, n_panels, 0, num_images, false);

    // Get shoeboxes for incident data
    for (std::size_t img_num = 0; img_num < num_images; ++img_num) {
      dxtbx::format::Image<double> img =
        incident_params.incident_data->get_corrected_data(img_num);
      dxtbx::format::Image<bool> mask =
        incident_params.incident_data->get_mask(img_num);

      dials::af::shared<scitbx::af::versa<double, scitbx::af::c_grid<2>>> output_data(
        n_panels);
      dials::af::shared<scitbx::af::versa<bool, scitbx::af::c_grid<2>>> output_mask(
        n_panels);

      for (std::size_t i = 0; i < output_data.size(); ++i) {
        output_data[i] = img.tile(i).data();
        output_mask[i] = mask.tile(i).data();
      }
      incident_shoebox_processor.next_data_only(
        dials::model::Image<double>(output_data.const_ref(), output_mask.const_ref()));
    }

    // Get shoeboxes for empty data
    for (std::size_t img_num = 0; img_num < num_images; ++img_num) {
      dxtbx::format::Image<double> img =
        incident_params.empty_data->get_corrected_data(img_num);
      dxtbx::format::Image<bool> mask = incident_params.empty_data->get_mask(img_num);

      dials::af::shared<scitbx::af::versa<double, scitbx::af::c_grid<2>>> output_data(
        n_panels);
      dials::af::shared<scitbx::af::versa<bool, scitbx::af::c_grid<2>>> output_mask(
        n_panels);

      for (std::size_t i = 0; i < output_data.size(); ++i) {
        output_data[i] = img.tile(i).data();
        output_mask[i] = mask.tile(i).data();
      }
      empty_shoebox_processor.next_data_only(
        dials::model::Image<double>(output_data.const_ref(), output_mask.const_ref()));
    }

    dials::af::shared<Shoebox<>> shoeboxes = reflection_table["shoebox"];
    dials::af::const_ref<int6> bboxes = reflection_table["bbox"];

    // Arrays to store ouput of integration
    dials::af::shared<bool> succeeded(reflection_table.size());
    dials::af::shared<double> intensities(reflection_table.size());
    dials::af::shared<double> variances(reflection_table.size());
    dials::af::shared<bool> succeeded_prf;
    dials::af::shared<double> intensities_prf;
    dials::af::shared<double> variances_prf;

    if (profile_params_1d || profile_params_3d) {
      succeeded_prf.resize(reflection_table.size());
      intensities_prf.resize(reflection_table.size());
      variances_prf.resize(reflection_table.size());
    }

    int bg_code = Valid | Background | BackgroundUsed;

    auto worker = [&](std::size_t start, std::size_t end) {
      for (std::size_t i = start; i < end; ++i) {
        Shoebox<> shoebox = shoeboxes[i];
        Shoebox<> i_shoebox = i_shoeboxes[i];
        Shoebox<> e_shoebox = e_shoeboxes[i];
        int6 bbox = bboxes[i];
        int panel = shoebox.panel;

        bool success = true;
        int n_background = 0;
        int n_signal = 0;
        double intensity = 0.0;
        double variance = 0.0;

        // First pass to get, n_signal, n_background, incident spectrum, empty spectrum
        // Shoebox data are ordered (z, y, x)
        scitbx::af::shared<double> incident_spectrum(shoebox.zsize(), 0.0);
        scitbx::af::shared<double> empty_spectrum(shoebox.zsize(), 0.0);
        std::vector<std::size_t> n_contrib(shoebox.zsize(), 0);

        for (std::size_t z = 0; z < shoebox.zsize(); ++z) {
          if (!success) {
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

              incident_spectrum[z] += i_shoebox.data(z, y, x);
              empty_spectrum[z] += e_shoebox.data(z, y, x);
              n_contrib[z]++;

              int mask = shoebox.mask(z, y, x);
              if ((mask & Foreground) == Foreground) {
                if ((mask & Valid) == Valid && (mask & Overlapped) == 0) {
                  n_signal++;
                } else {
                  success = false;
                }
              } else if ((mask & bg_code) == bg_code) {
                n_background++;
              }
            }
          }
        }

        // Smooth incident and empty to avoid dividing by a noisy signal
        scitbx::af::shared<double> smoothed_empty =
          dials_scaling::savitzky_golay(empty_spectrum, 7, 2);
        scitbx::af::shared<double> smoothed_incident =
          dials_scaling::savitzky_golay(incident_spectrum, 7, 2);

        // 1D profile fitting params
        scitbx::af::shared<double> tof_z(shoebox.zsize(), 0.0);
        scitbx::af::shared<double> projected_intensity(shoebox.zsize(), 0.0);

        // 3D profile fitting params done in xyz
        auto acc = shoebox.data.accessor();
        af::c_grid<3> transposed_acc(acc[2], acc[1], acc[0]);
        scitbx::af::versa<double, af::c_grid<3>> intensity_3d(transposed_acc);
        scitbx::af::versa<double, af::c_grid<3>> background_var_3d(transposed_acc);
        scitbx::af::versa<vec3<double>, af::c_grid<3>> coords_3d(transposed_acc);

        // Second pass to perform actual integration
        for (std::size_t z = 0; z < shoebox.zsize(); ++z) {
          double intensity_z = 0;

          if (!success) {
            break;
          }

          int frame_z = bbox[4] + z;
          double tof = img_tof[frame_z];
          tof_z[z] = tof;
          tof *= std::pow(10, -6);  // (s)

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

              // Raw counts
              double raw_S = shoebox.data(z, y, x);
              double raw_B = shoebox.background(z, y, x);
              double n = n_contrib[z];
              double raw_J = smoothed_incident[z] / n;
              double raw_E = smoothed_empty[z] / n;
              int mask = shoebox.mask(z, y, x);

              // Pixel data will be divided by incident spectrum
              if (raw_J < 1e-7) {
                continue;
              }

              // Normalize by proton charge
              double sample_pc = incident_params.sample_proton_charge;
              double incident_pc = incident_params.incident_proton_charge;
              double empty_pc = incident_params.empty_proton_charge;
              double S = raw_S / sample_pc;
              double B = raw_B / sample_pc;
              double J = raw_J / incident_pc;
              double E = raw_E / empty_pc;

              // Variances
              double var_S = std::abs(raw_S) / (sample_pc * sample_pc);
              double var_B = std::abs(raw_B) / (sample_pc * sample_pc);
              if (n_background > 0) {
                var_B *= (1.0 + double(n_signal) / double(n_background));
              }
              double var_J =
                (std::abs(raw_J) * n) / (incident_pc * incident_pc * n * n);
              double var_E = (std::abs(raw_E) * n) / (empty_pc * empty_pc * n * n);

              // Lorentz correction
              double two_theta = detector[panel].get_two_theta_at_pixel(
                unit_s0, scitbx::vec2<double>(panel_x, panel_y));
              scitbx::vec3<double> s1 = detector[panel].get_pixel_lab_coord(
                scitbx::vec2<double>(panel_x, panel_y));
              double distance = s1.length() + sample_to_source_distance;
              distance *= std::pow(10, -3);  // (m)
              double wl = ((Planck * tof) / (m_n * (distance))) * std::pow(10, 10);
              double sin_two_theta_sq = std::pow(sin(two_theta * .5), 2);
              double L = sin_two_theta_sq / std::pow(wl, 4);

              // Spherical absorption correction
              // for image data and incident data
              double two_theta_deg = two_theta * (180 / pi);
              int two_theta_idx = static_cast<int>(two_theta_deg / 10);
              double sample_muR =
                (corrections_data.sample_linear_scattering_c
                 + (corrections_data.sample_linear_absorption_c / 1.8) * wl)
                * corrections_data.sample_radius;
              double T = dials_scaling::tof_pixel_spherical_absorption_correction(
                sample_muR, two_theta, two_theta_idx);

              // Pixel data will be divided by the absorption correction
              if (T < 1e-7) {
                continue;
              }
              double incident_muR =
                (corrections_data.incident_linear_scattering_c
                 + (corrections_data.incident_linear_absorption_c / 1.8) * wl)
                * corrections_data.incident_radius;
              double J_T = dials_scaling::tof_pixel_spherical_absorption_correction(
                incident_muR, two_theta, two_theta_idx);

              // Pixel data will be divided by absorption correction
              if (J_T < 1e-7) {
                continue;
              }

              // Net signal
              double I0 = (S - B - E);
              double var_I0 = var_S + var_B + var_E;

              J /= J_T;

              if (J < 1e-7) {
                continue;
              }

              var_J /= (J_T * J_T);

              // Apply Lorentz correction and normalise by incident spectrum
              double I, var_I;
              if (apply_lorentz_correction) {
                I = L * I0 / (J * T);
                B *= L / (J * T);

                var_I = (L * L / (J * J * T * T)) * var_I0
                        + (L * L * I0 * I0 / (J * J * J * J * T * T)) * var_J;
                var_B = (L * L / (J * J * T * T)) * var_B
                        + (L * L * B * B / (J * J * J * J * T * T)) * var_J;
              } else {
                I = I0 / (J * T);
                B /= J * T;
                var_I = (1 / (J * J * T * T)) * var_I0
                        + (I0 * I0 / (J * J * J * J * T * T)) * var_J;
                var_B = (1 / (J * J * T * T)) * var_B
                        + (B * B / (J * J * J * J * T * T)) * var_J;
              }

              intensity_z += I;

              if (profile_params_3d) {
                intensity_3d(x, y, z) = I;
                background_var_3d(x, y, z) = var_B;
                double x_c = x + shoebox.xoffset() + 0.5;
                double y_c = y + shoebox.yoffset() + 0.5;
                double z_c = z + shoebox.zoffset() + 0.5;
                coords_3d(x, y, z) = vec3<double>(x_c, y_c, z_c);
              }

              // Accumulate if pixel in foreground & valid
              if ((mask & Foreground) == Foreground && (mask & Valid) == Valid
                  && (mask & Overlapped) == 0) {
                intensity += I;
                variance += var_I;
              } else if ((mask & Foreground) == Foreground) {
                success = false;
                break;
              }
            }
          }
          if (profile_params_1d) {
            projected_intensity[z] = intensity_z;
          }
        }

        // Overall values for shoebox summation
        succeeded[i] = success;
        intensities[i] = intensity;
        variances[i] = variance;

        if (profile_params_1d) {
          bool profile_success = false;

          double I_prf, var_prf;
          profile_success = fit_profile1d(projected_intensity.const_ref(),
                                          tof_z.const_ref(),
                                          *profile_params_1d,
                                          I_prf);
          if (profile_success) {
            intensities_prf[i] = I_prf;
            variances_prf[i] = variance;  // Use summation variance as approximation
          }
          succeeded_prf[i] = profile_success;
        } else if (profile_params_3d) {
          bool profile_success = false;

          double I_prf, var_prf;
          profile_success = fit_profile3d(coords_3d.const_ref(),
                                          intensity_3d,
                                          background_var_3d,
                                          *profile_params_3d,
                                          I_prf);
          if (profile_success) {
            intensities_prf[i] = I_prf;
            variances_prf[i] = variance;  // Use summation variance as approximation
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

    if (profile_params_1d || profile_params_3d) {
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
      if (profile_params_1d || profile_params_3d) {
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

  boost::python::tuple calculate_line_profile_for_reflection(
    dials::af::reflection_table &reflection,
    Experiment &experiment,
    ImageSequence &data,
    scitbx::af::shared<double> raw_projected_intensity_out,
    scitbx::af::shared<double> projected_intensity_out,
    scitbx::af::shared<double> projected_background_out,
    scitbx::af::shared<double> tof_z_out,
    const bool &apply_lorentz_correction) {
    /*
     * Calculates raw_projected_intensity, projected_intensity,
     * projected_background, sum_intensity, sum_variance
     */

    // This is a slight hack to make using current interfaces eaiser
    // E.g. the shoebox processor
    DIALS_ASSERT(reflection.size() == 1);

    Detector detector = *experiment.get_detector();
    Scan scan = *experiment.get_scan();

    // Required beam params
    std::shared_ptr<dxtbx::model::BeamBase> beam_ptr = experiment.get_beam();
    std::shared_ptr<PolychromaticBeam> beam =
      std::dynamic_pointer_cast<PolychromaticBeam>(beam_ptr);
    DIALS_ASSERT(beam != nullptr);
    vec3<double> unit_s0 = beam->get_unit_s0();
    double sample_to_source_distance = beam->get_sample_to_source_distance();

    // Required scan params
    scitbx::af::shared<double> img_tof = scan.get_property<double>("time_of_flight");

    // Required detector params
    int num_images = data.size();
    vec2<std::size_t> image_size = detector[0].get_image_size();
    DIALS_ASSERT(num_images == img_tof.size());

    // Mask codes
    int bg_code = Valid | Background | BackgroundUsed;

    dials::af::shared<Shoebox<>> shoeboxes = reflection["shoebox"];
    dials::model::Shoebox<> shoebox = shoeboxes[0];

    DIALS_ASSERT(raw_projected_intensity_out.size() == shoebox.zsize());
    DIALS_ASSERT(projected_intensity_out.size() == shoebox.zsize());
    DIALS_ASSERT(projected_background_out.size() == shoebox.zsize());
    DIALS_ASSERT(tof_z_out.size() == shoebox.zsize());

    int6 bbox = shoebox.bbox;
    int panel = shoebox.panel;

    bool success = true;
    int n_background = 0;
    int n_signal = 0;
    double intensity = 0.0;
    double variance = 0.0;

    // First pass to get, n_signal, n_background
    // Shoebox data are ordered (z, y, x)
    for (std::size_t z = 0; z < shoebox.zsize(); ++z) {
      if (!success) {
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

          int mask = shoebox.mask(z, y, x);
          if ((mask & Foreground) == Foreground) {
            if ((mask & Valid) == Valid && (mask & Overlapped) == 0) {
              n_signal++;
            } else {
              success = false;
            }
          } else if ((mask & bg_code) == bg_code) {
            n_background++;
          }
        }
      }
    }

    // Second pass to perform actual integration
    for (std::size_t z = 0; z < shoebox.zsize(); ++z) {
      double intensity_z = 0;
      double intensity_raw_z = 0;
      double background_z = 0;

      if (!success) {
        break;
      }

      int frame_z = bbox[4] + z;
      double tof = img_tof[frame_z];
      tof_z_out[z] = tof;
      tof *= std::pow(10, -6);  // (s)

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

          // Raw counts
          double raw_S = shoebox.data(z, y, x);
          double raw_B = shoebox.background(z, y, x);

          int mask = shoebox.mask(z, y, x);

          // Variances
          double var_S = std::abs(raw_S);
          double var_B = std::abs(raw_B);
          if (n_background > 0) {
            var_B *= (1.0 + double(n_signal) / double(n_background));
          }

          // Lorentz correction
          scitbx::vec3<double> s1 =
            detector[panel].get_pixel_lab_coord(scitbx::vec2<double>(panel_x, panel_y));
          double distance = s1.length() + sample_to_source_distance;
          distance *= std::pow(10, -3);  // (m)
          double wl = ((Planck * tof) / (m_n * (distance))) * std::pow(10, 10);
          double two_theta = detector[panel].get_two_theta_at_pixel(
            unit_s0, scitbx::vec2<double>(panel_x, panel_y));
          double sin_two_theta_sq = std::pow(sin(two_theta * .5), 2);
          double L = sin_two_theta_sq / std::pow(wl, 4);

          // Acount for some data having different sized bin widths
          double bin_width_correction = img_tof[frame_z] - img_tof[frame_z - 1];

          // Net signal
          double I0 = raw_S - raw_B;
          double var_I0 = var_S + var_B;
          double I_raw = raw_S;

          double I = I0 / bin_width_correction;
          double B = raw_B / bin_width_correction;
          I_raw /= bin_width_correction;
          double var_I = var_I0 / (bin_width_correction * bin_width_correction);

          if (apply_lorentz_correction) {
            I *= L;
            var_I *= (L * L);
            I_raw *= L;
            B *= L;
          }

          intensity_z += I;
          intensity_raw_z += I_raw;
          background_z += B;

          // Accumulate if pixel in foreground & valid
          if ((mask & Foreground) == Foreground && (mask & Valid) == Valid
              && (mask & Overlapped) == 0) {
            intensity += I;
            variance += var_I;
          } else if ((mask & Foreground) == Foreground) {
            success = false;
            break;
          }
        }
      }
      raw_projected_intensity_out[z] = intensity_raw_z;
      projected_intensity_out[z] = intensity_z;
      projected_background_out[z] = background_z;
    }

    return boost::python::make_tuple(intensity, variance, success);
  }

  boost::python::tuple calculate_line_profile_for_reflection_3d(
    dials::af::reflection_table &reflection,
    Experiment &experiment,
    ImageSequence &data,
    scitbx::af::shared<vec3<double>> coords,
    scitbx::af::shared<double> raw_projected_intensity_out,
    scitbx::af::shared<double> projected_intensity_out,
    scitbx::af::shared<double> projected_background_out,
    scitbx::af::shared<double> tof_z_out,
    const bool &apply_lorentz_correction,
    TOFProfile3DParams &profile_params_3d) {
    /*
     * Calculates raw_projected_intensity, projected_intensity,
     * projected_background, sum_intensity, sum_variance
     */

    // This is a slight hack to make using current interfaces eaiser
    // E.g. the shoebox processor
    DIALS_ASSERT(reflection.size() == 1);

    Detector detector = *experiment.get_detector();
    Scan scan = *experiment.get_scan();

    // Required beam params
    std::shared_ptr<dxtbx::model::BeamBase> beam_ptr = experiment.get_beam();
    std::shared_ptr<PolychromaticBeam> beam =
      std::dynamic_pointer_cast<PolychromaticBeam>(beam_ptr);
    DIALS_ASSERT(beam != nullptr);
    vec3<double> unit_s0 = beam->get_unit_s0();
    double sample_to_source_distance = beam->get_sample_to_source_distance();

    // Required scan params
    scitbx::af::shared<double> img_tof = scan.get_property<double>("time_of_flight");

    // Required detector params
    int num_images = data.size();
    vec2<std::size_t> image_size = detector[0].get_image_size();
    DIALS_ASSERT(num_images == img_tof.size());

    // Mask codes
    int bg_code = Valid | Background | BackgroundUsed;

    dials::af::shared<Shoebox<>> shoeboxes = reflection["shoebox"];
    dials::model::Shoebox<> shoebox = shoeboxes[0];

    DIALS_ASSERT(raw_projected_intensity_out.size() == shoebox.zsize());
    DIALS_ASSERT(projected_intensity_out.size() == shoebox.zsize());
    DIALS_ASSERT(projected_background_out.size() == shoebox.zsize());
    DIALS_ASSERT(tof_z_out.size() == shoebox.zsize());

    int6 bbox = shoebox.bbox;
    int panel = shoebox.panel;

    bool sum_success = true;
    int n_background = 0;
    int n_signal = 0;
    double I_sum = 0.0;
    double var_sum = 0.0;

    // First pass to get, n_signal, n_background
    // Shoebox data are ordered (z, y, x)
    for (std::size_t z = 0; z < shoebox.zsize(); ++z) {
      if (!sum_success) {
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

          int mask = shoebox.mask(z, y, x);
          if ((mask & Foreground) == Foreground) {
            if ((mask & Valid) == Valid && (mask & Overlapped) == 0) {
              n_signal++;
            } else {
              sum_success = false;
            }
          } else if ((mask & bg_code) == bg_code) {
            n_background++;
          }
        }
      }
    }
    // 3D profile fitting params done in xyz
    auto acc = shoebox.data.accessor();
    af::c_grid<3> transposed_acc(acc[2], acc[1], acc[0]);
    scitbx::af::versa<double, af::c_grid<3>> intensity_3d(transposed_acc);
    scitbx::af::versa<double, af::c_grid<3>> background_var_3d(transposed_acc);
    scitbx::af::versa<vec3<double>, af::c_grid<3>> coords_3d(transposed_acc);
    int coord_count = 0;

    // Second pass to perform actual integration
    for (std::size_t z = 0; z < shoebox.zsize(); ++z) {
      double intensity_z = 0;
      double intensity_raw_z = 0;
      double background_z = 0;

      if (!sum_success) {
        break;
      }

      int frame_z = bbox[4] + z;
      double tof = img_tof[frame_z];
      tof_z_out[z] = tof;
      tof *= std::pow(10, -6);  // (s)

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

          // Raw counts
          double raw_S = shoebox.data(z, y, x);
          double raw_B = shoebox.background(z, y, x);

          int mask = shoebox.mask(z, y, x);

          // Variances
          double var_S = std::abs(raw_S);
          double var_B = std::abs(raw_B);
          if (n_background > 0) {
            var_B *= (1.0 + double(n_signal) / double(n_background));
          }

          // Lorentz correction
          scitbx::vec3<double> s1 =
            detector[panel].get_pixel_lab_coord(scitbx::vec2<double>(panel_x, panel_y));
          double distance = s1.length() + sample_to_source_distance;
          distance *= std::pow(10, -3);  // (m)
          double wl = ((Planck * tof) / (m_n * (distance))) * std::pow(10, 10);
          double two_theta = detector[panel].get_two_theta_at_pixel(
            unit_s0, scitbx::vec2<double>(panel_x, panel_y));
          double sin_two_theta_sq = std::pow(sin(two_theta * .5), 2);
          double L = sin_two_theta_sq / std::pow(wl, 4);

          // Acount for some data having different sized bin widths
          double bin_width_correction = img_tof[frame_z] - img_tof[frame_z - 1];

          // Net signal
          double I0 = raw_S - raw_B;
          double var_I0 = var_S + var_B;
          double I_raw = raw_S;
          double B = raw_B;

          double I = I0 / bin_width_correction;
          I_raw /= bin_width_correction;
          double var_I = var_I0 / (bin_width_correction * bin_width_correction);
          B /= bin_width_correction;
          var_B /= (bin_width_correction * bin_width_correction);

          if (apply_lorentz_correction) {
            I *= L;
            var_I *= (L * L);
            var_B *= (L * L);
            I_raw *= L;
            B *= L;
          }

          intensity_z += I;
          intensity_raw_z += I_raw;
          background_z += B;

          intensity_3d(x, y, z) = I;
          background_var_3d(x, y, z) = var_B;
          coords_3d(x, y, z) = coords[coord_count];
          coord_count++;

          // Accumulate if pixel in foreground & valid
          if ((mask & Foreground) == Foreground && (mask & Valid) == Valid
              && (mask & Overlapped) == 0) {
            I_sum += I;
            var_sum += var_I;
          } else if ((mask & Foreground) == Foreground) {
            sum_success = false;
            break;
          }
        }
      }
      raw_projected_intensity_out[z] = intensity_raw_z;
      projected_intensity_out[z] = intensity_z;
      projected_background_out[z] = background_z;
    }

    double I_prf = 0;
    double var_prf = 0;
    bool profile_success = false;
    scitbx::af::versa<double, scitbx::af::c_grid<3>> profile_3d_out(
      intensity_3d.accessor());

    if (sum_success) {
      profile_success = fit_profile3d(coords_3d.const_ref(),
                                      intensity_3d,
                                      background_var_3d,
                                      profile_params_3d,
                                      I_prf,
                                      profile_3d_out);
    } else {
      profile_success = false;
    }
    return boost::python::make_tuple(
      I_prf, var_prf, I_sum, var_sum, profile_success, profile_3d_out);
  }

  boost::python::tuple calculate_line_profile_for_reflection(
    dials::af::reflection_table &reflection,
    Experiment &experiment,
    ImageSequence &data,
    scitbx::af::shared<double> raw_projected_intensity_out,
    scitbx::af::shared<double> projected_intensity_out,
    scitbx::af::shared<double> projected_background_out,
    scitbx::af::shared<double> tof_z_out,
    scitbx::af::shared<double> line_profile_out,
    const bool &apply_lorentz_correction,
    TOFProfile1DParams &profile_params_1d) {
    /*
     * Calculates raw_projected_intensity, projected_intensity, line_profile
     * projected_background, sum_intensity, sum_variance, prf_intensity, prf_variance
     */

    boost::python::tuple result =
      calculate_line_profile_for_reflection(reflection,
                                            experiment,
                                            data,
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
      profile_success = fit_profile1d(projected_intensity_out.const_ref(),
                                      tof_z_out.const_ref(),
                                      profile_params_1d,
                                      I_prf,
                                      line_profile_out);

    } else {
      profile_success = false;
    }
    return boost::python::make_tuple(I_prf, var_prf, I_sum, var_sum, profile_success);
  }

  boost::python::tuple calculate_line_profile_for_reflection(
    dials::af::reflection_table &reflection,
    Experiment &experiment,
    ImageSequence &data,
    const dials_scaling::TOFIncidentSpectrumParams &incident_params,
    scitbx::af::shared<double> raw_projected_intensity_out,
    scitbx::af::shared<double> projected_intensity_out,
    scitbx::af::shared<double> projected_background_out,
    scitbx::af::shared<double> tof_z_out,
    const bool &apply_lorentz_correction) {
    /*
     * Calculates raw_projected_intensity, projected_intensity,
     * projected_background, sum_intensity, sum_variance
     */

    // This is a slight hack to make using current interfaces eaiser
    // E.g. the shoebox processor
    DIALS_ASSERT(reflection.size() == 1);
    dials::af::shared<Shoebox<>> shoeboxes = reflection["shoebox"];
    dials::model::Shoebox<> shoebox = shoeboxes[0];

    Detector detector = *experiment.get_detector();
    Scan scan = *experiment.get_scan();

    // Required beam params
    std::shared_ptr<dxtbx::model::BeamBase> beam_ptr = experiment.get_beam();
    std::shared_ptr<PolychromaticBeam> beam =
      std::dynamic_pointer_cast<PolychromaticBeam>(beam_ptr);
    DIALS_ASSERT(beam != nullptr);
    vec3<double> unit_s0 = beam->get_unit_s0();
    double sample_to_source_distance = beam->get_sample_to_source_distance();

    // Required scan params
    scitbx::af::shared<double> img_tof = scan.get_property<double>("time_of_flight");

    // Required detector params
    int num_images = data.size();
    vec2<std::size_t> image_size = detector[0].get_image_size();
    DIALS_ASSERT(num_images == img_tof.size());
    int n_panels = detector.size();

    int6 bbox = shoebox.bbox;
    int panel = shoebox.panel;

    // Copy reflections for incident and empty runs
    boost::python::dict d;
    dials::af::reflection_table i_reflection_table =
      dxtbx::af::flex_table_suite::deepcopy(reflection, d);
    dials::af::reflection_table e_reflection_table =
      dxtbx::af::flex_table_suite::deepcopy(reflection, d);

    dials::af::ref<Shoebox<>> i_shoeboxes = i_reflection_table["shoebox"];
    dials::af::ref<Shoebox<>> e_shoeboxes = e_reflection_table["shoebox"];
    i_shoeboxes[0].deallocate();
    e_shoeboxes[0].deallocate();

    ShoeboxProcessor incident_shoebox_processor(
      i_reflection_table, n_panels, 0, num_images, false);

    ShoeboxProcessor empty_shoebox_processor(
      e_reflection_table, n_panels, 0, num_images, false);

    // Get shoebox for incident data
    for (std::size_t img_num = 0; img_num < num_images; ++img_num) {
      dxtbx::format::Image<double> img =
        incident_params.incident_data->get_corrected_data(img_num);
      dxtbx::format::Image<bool> mask =
        incident_params.incident_data->get_mask(img_num);

      dials::af::shared<scitbx::af::versa<double, scitbx::af::c_grid<2>>> output_data(
        n_panels);
      dials::af::shared<scitbx::af::versa<bool, scitbx::af::c_grid<2>>> output_mask(
        n_panels);

      for (std::size_t i = 0; i < output_data.size(); ++i) {
        output_data[i] = img.tile(i).data();
        output_mask[i] = mask.tile(i).data();
      }
      incident_shoebox_processor.next_data_only(
        dials::model::Image<double>(output_data.const_ref(), output_mask.const_ref()));
    }

    // Get shoeboxes for empty data
    for (std::size_t img_num = 0; img_num < num_images; ++img_num) {
      dxtbx::format::Image<double> img =
        incident_params.empty_data->get_corrected_data(img_num);
      dxtbx::format::Image<bool> mask = incident_params.empty_data->get_mask(img_num);

      dials::af::shared<scitbx::af::versa<double, scitbx::af::c_grid<2>>> output_data(
        n_panels);
      dials::af::shared<scitbx::af::versa<bool, scitbx::af::c_grid<2>>> output_mask(
        n_panels);

      for (std::size_t i = 0; i < output_data.size(); ++i) {
        output_data[i] = img.tile(i).data();
        output_mask[i] = mask.tile(i).data();
      }
      empty_shoebox_processor.next_data_only(
        dials::model::Image<double>(output_data.const_ref(), output_mask.const_ref()));
    }

    // Mask codes
    int bg_code = Valid | Background | BackgroundUsed;

    bool success = true;
    int n_background = 0;
    int n_signal = 0;
    double intensity = 0.0;
    double variance = 0.0;

    scitbx::af::shared<double> incident_spectrum(shoebox.zsize(), 0.0);
    scitbx::af::shared<double> empty_spectrum(shoebox.zsize(), 0.0);
    std::vector<std::size_t> n_contrib(shoebox.zsize(), 0);
    Shoebox<> i_shoebox = i_shoeboxes[0];
    Shoebox<> e_shoebox = e_shoeboxes[0];

    // First pass to get, n_signal, n_background
    // Shoebox data are ordered (z, y, x)
    for (std::size_t z = 0; z < shoebox.zsize(); ++z) {
      if (!success) {
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

          incident_spectrum[z] += i_shoebox.data(z, y, x);
          empty_spectrum[z] += e_shoebox.data(z, y, x);
          n_contrib[z]++;

          int mask = shoebox.mask(z, y, x);
          if ((mask & Foreground) == Foreground) {
            if ((mask & Valid) == Valid && (mask & Overlapped) == 0) {
              n_signal++;
            } else {
              success = false;
            }
          } else if ((mask & bg_code) == bg_code) {
            n_background++;
          }
        }
      }
    }
    // Smooth incident and empty to avoid dividing by a noisy signal
    scitbx::af::shared<double> smoothed_empty =
      dials_scaling::savitzky_golay(empty_spectrum, 7, 2);
    scitbx::af::shared<double> smoothed_incident =
      dials_scaling::savitzky_golay(incident_spectrum, 7, 2);

    DIALS_ASSERT(raw_projected_intensity_out.size() == shoebox.zsize());
    DIALS_ASSERT(projected_intensity_out.size() == shoebox.zsize());
    DIALS_ASSERT(projected_background_out.size() == shoebox.zsize());
    DIALS_ASSERT(tof_z_out.size() == shoebox.zsize());

    // Second pass to perform actual integration
    for (std::size_t z = 0; z < shoebox.zsize(); ++z) {
      double intensity_z = 0;
      double intensity_raw_z = 0;
      double background_z = 0;

      if (!success) {
        break;
      }

      int frame_z = bbox[4] + z;
      double tof = img_tof[frame_z];
      tof_z_out[z] = tof;
      tof *= std::pow(10, -6);  // (s)

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
          // Raw counts
          double raw_S = shoebox.data(z, y, x);
          double raw_B = shoebox.background(z, y, x);
          double n = n_contrib[z];
          double raw_J = smoothed_incident[z] / n;
          double raw_E = smoothed_empty[z] / n;

          int mask = shoebox.mask(z, y, x);

          // Normalise by proton charge
          double sample_pc = incident_params.sample_proton_charge;
          double incident_pc = incident_params.incident_proton_charge;
          double empty_pc = incident_params.empty_proton_charge;
          double S = raw_S / sample_pc;
          double B = raw_B / sample_pc;
          double J = raw_J / incident_pc;
          double E = raw_E / empty_pc;

          // Variances
          double var_S = std::abs(raw_S) / (sample_pc * sample_pc);
          double var_B = std::abs(raw_B) / (sample_pc * sample_pc);
          if (n_background > 0) {
            var_B *= (1.0 + double(n_signal) / double(n_background));
          }
          double var_J = (std::abs(raw_J) * n) / (incident_pc * incident_pc * n * n);
          double var_E = (std::abs(raw_E) * n) / (empty_pc * empty_pc * n * n);

          // Lorentz correction
          scitbx::vec3<double> s1 =
            detector[panel].get_pixel_lab_coord(scitbx::vec2<double>(panel_x, panel_y));
          double distance = s1.length() + sample_to_source_distance;
          distance *= std::pow(10, -3);  // (m)
          double wl = ((Planck * tof) / (m_n * (distance))) * std::pow(10, 10);
          double two_theta = detector[panel].get_two_theta_at_pixel(
            unit_s0, scitbx::vec2<double>(panel_x, panel_y));
          double sin_two_theta_sq = std::pow(sin(two_theta * .5), 2);
          double L = sin_two_theta_sq / std::pow(wl, 4);

          // Net signal
          double I0 = S - B - E;
          double var_I0 = var_S + var_B + var_E;
          double I_raw = S;

          // Apply Lorentz correction and normalise by incident spectrum
          double I, var_I;
          if (apply_lorentz_correction) {
            I = L * I0 / J;
            I_raw *= L / J;
            B *= L / J;
            // Var( L*I0/J ) = (L/J)^2 Var(I0) + (L*I0/J^2)^2 Var(J)
            var_I =
              (L * L / (J * J)) * var_I0 + (L * L * I0 * I0 / (J * J * J * J)) * var_J;
          } else {
            I = I0 / J;
            I_raw /= J;
            B /= J;
            // Var( I0/J ) = (1/J)^2 Var(I0) + (I0/J^2)^2 Var(J)
            var_I = (1 / (J * J)) * var_I0 + (I0 * I0 / (J * J * J * J)) * var_J;
          }

          intensity_z += I;
          intensity_raw_z += I_raw;
          background_z += B;

          // Accumulate if pixel in foreground & valid
          if ((mask & Foreground) == Foreground && (mask & Valid) == Valid
              && (mask & Overlapped) == 0) {
            intensity += I;
            variance += var_I;
          } else if ((mask & Foreground) == Foreground) {
            success = false;
            break;
          }
        }
      }
      raw_projected_intensity_out[z] = intensity_raw_z;
      projected_intensity_out[z] = intensity_z;
      projected_background_out[z] = background_z;
    }

    return boost::python::make_tuple(intensity, variance, success);
  }

  boost::python::tuple calculate_line_profile_for_reflection(
    dials::af::reflection_table &reflection,
    Experiment &experiment,
    ImageSequence &data,
    const dials_scaling::TOFIncidentSpectrumParams &incident_params,
    scitbx::af::shared<double> raw_projected_intensity_out,
    scitbx::af::shared<double> projected_intensity_out,
    scitbx::af::shared<double> projected_background_out,
    scitbx::af::shared<double> tof_z_out,
    scitbx::af::shared<double> line_profile_out,
    const bool &apply_lorentz_correction,
    TOFProfile1DParams &profile_params_1d) {
    /*
     * Calculates raw_projected_intensity, projected_intensity, line_profile
     * projected_background, sum_intensity, sum_variance, prf_intensity, prf_variance
     */

    boost::python::tuple result =
      calculate_line_profile_for_reflection(reflection,
                                            experiment,
                                            data,
                                            incident_params,
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
      profile_success = fit_profile1d(projected_intensity_out.const_ref(),
                                      tof_z_out.const_ref(),
                                      profile_params_1d,
                                      I_prf,
                                      line_profile_out);

    } else {
      profile_success = false;
    }
    return boost::python::make_tuple(I_prf, var_prf, I_sum, var_sum, profile_success);
  }

  boost::python::tuple calculate_line_profile_for_reflection(
    dials::af::reflection_table &reflection,
    Experiment &experiment,
    ImageSequence &data,
    const dials_scaling::TOFIncidentSpectrumParams &incident_params,
    const dials_scaling::TOFAbsorptionParams &corrections_data,
    scitbx::af::shared<double> raw_projected_intensity_out,
    scitbx::af::shared<double> projected_intensity_out,
    scitbx::af::shared<double> projected_background_out,
    scitbx::af::shared<double> tof_z_out,
    const bool &apply_lorentz_correction) {
    /*
     * Calculates raw_projected_intensity, projected_intensity,
     * projected_background, sum_intensity, sum_variance
     */

    // This is a slight hack to make using current interfaces eaiser
    // E.g. the shoebox processor
    DIALS_ASSERT(reflection.size() == 1);
    dials::af::shared<Shoebox<>> shoeboxes = reflection["shoebox"];
    dials::model::Shoebox<> shoebox = shoeboxes[0];

    Detector detector = *experiment.get_detector();
    Scan scan = *experiment.get_scan();

    // Required beam params
    std::shared_ptr<dxtbx::model::BeamBase> beam_ptr = experiment.get_beam();
    std::shared_ptr<PolychromaticBeam> beam =
      std::dynamic_pointer_cast<PolychromaticBeam>(beam_ptr);
    DIALS_ASSERT(beam != nullptr);
    vec3<double> unit_s0 = beam->get_unit_s0();
    double sample_to_source_distance = beam->get_sample_to_source_distance();

    // Required scan params
    scitbx::af::shared<double> img_tof = scan.get_property<double>("time_of_flight");

    // Required detector params
    int num_images = data.size();
    vec2<std::size_t> image_size = detector[0].get_image_size();
    DIALS_ASSERT(num_images == img_tof.size());
    int n_panels = detector.size();

    int6 bbox = shoebox.bbox;
    int panel = shoebox.panel;

    // Copy reflections for incident and empty runs
    boost::python::dict d;
    dials::af::reflection_table i_reflection_table =
      dxtbx::af::flex_table_suite::deepcopy(reflection, d);
    dials::af::reflection_table e_reflection_table =
      dxtbx::af::flex_table_suite::deepcopy(reflection, d);

    dials::af::ref<Shoebox<>> i_shoeboxes = i_reflection_table["shoebox"];
    dials::af::ref<Shoebox<>> e_shoeboxes = e_reflection_table["shoebox"];
    i_shoeboxes[0].deallocate();
    e_shoeboxes[0].deallocate();

    ShoeboxProcessor incident_shoebox_processor(
      i_reflection_table, n_panels, 0, num_images, false);

    ShoeboxProcessor empty_shoebox_processor(
      e_reflection_table, n_panels, 0, num_images, false);

    // Get shoebox for incident data
    for (std::size_t img_num = 0; img_num < num_images; ++img_num) {
      dxtbx::format::Image<double> img =
        incident_params.incident_data->get_corrected_data(img_num);
      dxtbx::format::Image<bool> mask =
        incident_params.incident_data->get_mask(img_num);

      dials::af::shared<scitbx::af::versa<double, scitbx::af::c_grid<2>>> output_data(
        n_panels);
      dials::af::shared<scitbx::af::versa<bool, scitbx::af::c_grid<2>>> output_mask(
        n_panels);

      for (std::size_t i = 0; i < output_data.size(); ++i) {
        output_data[i] = img.tile(i).data();
        output_mask[i] = mask.tile(i).data();
      }
      incident_shoebox_processor.next_data_only(
        dials::model::Image<double>(output_data.const_ref(), output_mask.const_ref()));
    }

    // Get shoeboxes for empty data
    for (std::size_t img_num = 0; img_num < num_images; ++img_num) {
      dxtbx::format::Image<double> img =
        incident_params.empty_data->get_corrected_data(img_num);
      dxtbx::format::Image<bool> mask = incident_params.empty_data->get_mask(img_num);

      dials::af::shared<scitbx::af::versa<double, scitbx::af::c_grid<2>>> output_data(
        n_panels);
      dials::af::shared<scitbx::af::versa<bool, scitbx::af::c_grid<2>>> output_mask(
        n_panels);

      for (std::size_t i = 0; i < output_data.size(); ++i) {
        output_data[i] = img.tile(i).data();
        output_mask[i] = mask.tile(i).data();
      }
      empty_shoebox_processor.next_data_only(
        dials::model::Image<double>(output_data.const_ref(), output_mask.const_ref()));
    }

    // Mask codes
    int bg_code = Valid | Background | BackgroundUsed;

    bool success = true;
    int n_background = 0;
    int n_signal = 0;
    double intensity = 0.0;
    double variance = 0.0;

    scitbx::af::shared<double> incident_spectrum(shoebox.zsize(), 0.0);
    scitbx::af::shared<double> empty_spectrum(shoebox.zsize(), 0.0);
    std::vector<std::size_t> n_contrib(shoebox.zsize(), 0);
    dials::model::Shoebox<> i_shoebox = i_shoeboxes[0];
    dials::model::Shoebox<> e_shoebox = e_shoeboxes[0];

    // First pass to get, n_signal, n_background
    // Shoebox data are ordered (z, y, x)
    for (std::size_t z = 0; z < shoebox.zsize(); ++z) {
      if (!success) {
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

          incident_spectrum[z] += i_shoebox.data(z, y, x);
          empty_spectrum[z] += e_shoebox.data(z, y, x);
          n_contrib[z]++;

          int mask = shoebox.mask(z, y, x);
          if ((mask & Foreground) == Foreground) {
            if ((mask & Valid) == Valid && (mask & Overlapped) == 0) {
              n_signal++;
            } else {
              success = false;
            }
          } else if ((mask & bg_code) == bg_code) {
            n_background++;
          }
        }
      }
    }

    DIALS_ASSERT(raw_projected_intensity_out.size() == shoebox.zsize());
    DIALS_ASSERT(projected_intensity_out.size() == shoebox.zsize());
    DIALS_ASSERT(projected_background_out.size() == shoebox.zsize());
    DIALS_ASSERT(tof_z_out.size() == shoebox.zsize());

    // Smooth incident and empty to avoid dividing by a noisy signal
    scitbx::af::shared<double> smoothed_empty =
      dials_scaling::savitzky_golay(empty_spectrum, 7, 2);
    scitbx::af::shared<double> smoothed_incident =
      dials_scaling::savitzky_golay(incident_spectrum, 7, 2);

    // Second pass to perform actual integration
    for (std::size_t z = 0; z < shoebox.zsize(); ++z) {
      double intensity_z = 0;
      double intensity_raw_z = 0;
      double background_z = 0;

      if (!success) {
        break;
      }

      int frame_z = bbox[4] + z;
      double tof = img_tof[frame_z];
      tof_z_out[z] = tof;
      tof *= std::pow(10, -6);  // (s)

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
          // Raw counts
          double raw_S = shoebox.data(z, y, x);
          double raw_B = shoebox.background(z, y, x);
          double n = n_contrib[z];
          double raw_J = smoothed_incident[z] / n;
          double raw_E = smoothed_empty[z] / n;

          int mask = shoebox.mask(z, y, x);

          // Normalise by proton charge
          double sample_pc = incident_params.sample_proton_charge;
          double incident_pc = incident_params.incident_proton_charge;
          double empty_pc = incident_params.empty_proton_charge;
          double S = raw_S / sample_pc;
          double B = raw_B / sample_pc;
          double J = raw_J / incident_pc;
          double E = raw_E / empty_pc;

          // Variances
          double var_S = std::abs(raw_S) / (sample_pc * sample_pc);
          double var_B = std::abs(raw_B) / (sample_pc * sample_pc);
          if (n_background > 0) {
            var_B *= (1.0 + double(n_signal) / double(n_background));
          }
          double var_J = (std::abs(raw_J) * n) / (incident_pc * incident_pc * n * n);
          double var_E = (std::abs(raw_E) * n) / (empty_pc * empty_pc * n * n);

          // Lorentz correction
          scitbx::vec3<double> s1 =
            detector[panel].get_pixel_lab_coord(scitbx::vec2<double>(panel_x, panel_y));
          double distance = s1.length() + sample_to_source_distance;
          distance *= std::pow(10, -3);  // (m)
          double wl = ((Planck * tof) / (m_n * (distance))) * std::pow(10, 10);
          double two_theta = detector[panel].get_two_theta_at_pixel(
            unit_s0, scitbx::vec2<double>(panel_x, panel_y));
          double sin_two_theta_sq = std::pow(sin(two_theta * .5), 2);
          double L = sin_two_theta_sq / std::pow(wl, 4);

          // Spherical absorption correction
          // for image data and incident data
          double two_theta_deg = two_theta * (180 / pi);
          int two_theta_idx = static_cast<int>(two_theta_deg / 10);
          double sample_muR =
            (corrections_data.sample_linear_scattering_c
             + (corrections_data.sample_linear_absorption_c / 1.8) * wl)
            * corrections_data.sample_radius;
          double T = dials_scaling::tof_pixel_spherical_absorption_correction(
            sample_muR, two_theta, two_theta_idx);

          // Pixel data will be divided by the absorption correction
          if (T < 1e-7) {
            continue;
          }
          double incident_muR =
            (corrections_data.incident_linear_scattering_c
             + (corrections_data.incident_linear_absorption_c / 1.8) * wl)
            * corrections_data.incident_radius;
          double J_T = dials_scaling::tof_pixel_spherical_absorption_correction(
            incident_muR, two_theta, two_theta_idx);

          // Pixel data will be divided by absorption correction
          if (J_T < 1e-7) {
            continue;
          }

          // Net signal
          double I0 = S - B - E;
          double var_I0 = var_S + var_B + var_E;
          double I_raw = raw_S;

          J /= J_T;

          if (J < 1e-7) {
            continue;
          }

          var_J /= (J_T * J_T);

          // Apply Lorentz correction and normalise by incident spectrum
          double I, var_I;
          if (apply_lorentz_correction) {
            I = L * I0 / (J * T);
            I_raw *= (L / (J * T));
            B *= (L / (J * T));

            var_I = (L * L / (J * J * T * T)) * var_I0
                    + (L * L * I0 * I0 / (J * J * J * J * T * T)) * var_J;
          } else {
            I = I0 / (J * T);
            var_I = (1 / (J * J * T * T)) * var_I0
                    + (I0 * I0 / (J * J * J * J * T * T)) * var_J;
            I_raw /= (J * T);
            B /= (J * T);
          }

          intensity_z += I;
          intensity_raw_z += I_raw;
          background_z += B;

          // Accumulate if pixel in foreground & valid
          if ((mask & Foreground) == Foreground && (mask & Valid) == Valid
              && (mask & Overlapped) == 0) {
            intensity += I;
            variance += var_I;
          } else if ((mask & Foreground) == Foreground) {
            success = false;
            break;
          }
        }
      }
      raw_projected_intensity_out[z] = intensity_raw_z;
      projected_intensity_out[z] = intensity_z;
      projected_background_out[z] = background_z;
    }

    return boost::python::make_tuple(intensity, variance, success);
  }

  boost::python::tuple calculate_line_profile_for_reflection(
    dials::af::reflection_table &reflection,
    Experiment &experiment,
    ImageSequence &data,
    const dials_scaling::TOFIncidentSpectrumParams &incident_params,
    const dials_scaling::TOFAbsorptionParams &corrections_data,
    scitbx::af::shared<double> raw_projected_intensity_out,
    scitbx::af::shared<double> projected_intensity_out,
    scitbx::af::shared<double> projected_background_out,
    scitbx::af::shared<double> tof_z_out,
    scitbx::af::shared<double> line_profile_out,
    const bool &apply_lorentz_correction,
    TOFProfile1DParams &profile_params_1d) {
    /*
     * Calculates raw_projected_intensity, projected_intensity, line_profile
     * projected_background, sum_intensity, sum_variance, prf_intensity, prf_variance
     */

    boost::python::tuple result =
      calculate_line_profile_for_reflection(reflection,
                                            experiment,
                                            data,
                                            incident_params,
                                            corrections_data,
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
      profile_success = fit_profile1d(projected_intensity_out.const_ref(),
                                      tof_z_out.const_ref(),
                                      profile_params_1d,
                                      I_prf,
                                      line_profile_out);

    } else {
      profile_success = false;
    }
    return boost::python::make_tuple(I_prf, var_prf, I_sum, var_sum, profile_success);
  }

}}  // namespace dials::algorithms

#endif /* DIALS_ALGORITHMS_INTEGRATION_TOF_INTEGRATION_H */