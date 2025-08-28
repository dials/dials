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

  // Taken from
  // https://github.com/mantidproject/mantid/blob/30e97650e69dec5f6edbc456aa81f2d8f1715fa3/Framework/Crystal/inc/MantidCrystal/AnvredCorrection.h
  // Based on the method outlined in
  // C.W Dwiggins Jnr, Rapid calculation of X-ray absorption correction factors for
  // spheres to an accuracy of 0.05%, Acta Cryst. 1975, A31, 395-396
  // https://doi.org/10.1107/S0567739475000873
  const double pc[8][19] = {{-6.4910e-07,
                             -6.8938e-07,
                             -7.8149e-07,
                             8.1682e-08,
                             1.8008e-06,
                             3.3916e-06,
                             4.5095e-06,
                             4.7970e-06,
                             4.4934e-06,
                             3.6700e-06,
                             2.5881e-06,
                             1.5007e-06,
                             3.7669e-07,
                             -7.9487e-07,
                             -1.7935e-06,
                             -2.5563e-06,
                             -3.1113e-06,
                             -3.3993e-06,
                             -3.5091e-06},
                            {1.0839e-05,
                             1.1582e-05,
                             1.1004e-05,
                             -2.2848e-05,
                             -8.1974e-05,
                             -1.3268e-04,
                             -1.6486e-04,
                             -1.6839e-04,
                             -1.5242e-04,
                             -1.1949e-04,
                             -7.8682e-05,
                             -3.7973e-05,
                             2.9117e-06,
                             4.4823e-05,
                             8.0464e-05,
                             1.0769e-04,
                             1.2753e-04,
                             1.3800e-04,
                             1.4190e-04},
                            {8.7140e-05,
                             9.0870e-05,
                             1.6706e-04,
                             6.9008e-04,
                             1.4781e-03,
                             2.0818e-03,
                             2.3973e-03,
                             2.3209e-03,
                             1.9935e-03,
                             1.4508e-03,
                             8.1903e-04,
                             1.9608e-04,
                             -4.1128e-04,
                             -1.0205e-03,
                             -1.5374e-03,
                             -1.9329e-03,
                             -2.2212e-03,
                             -2.3760e-03,
                             -2.4324e-03},
                            {-2.9549e-03,
                             -3.1360e-03,
                             -4.2431e-03,
                             -8.1103e-03,
                             -1.2989e-02,
                             -1.6012e-02,
                             -1.6815e-02,
                             -1.4962e-02,
                             -1.1563e-02,
                             -6.8581e-03,
                             -1.7302e-03,
                             3.2400e-03,
                             7.9409e-03,
                             1.2528e-02,
                             1.6414e-02,
                             1.9394e-02,
                             2.1568e-02,
                             2.2758e-02,
                             2.3182e-02},
                            {1.7934e-02,
                             1.9304e-02,
                             2.4706e-02,
                             3.6759e-02,
                             4.8351e-02,
                             5.1049e-02,
                             4.5368e-02,
                             3.0864e-02,
                             1.2086e-02,
                             -1.0254e-02,
                             -3.2992e-02,
                             -5.4495e-02,
                             -7.4205e-02,
                             -9.2818e-02,
                             -1.0855e-01,
                             -1.2068e-01,
                             -1.2954e-01,
                             -1.3451e-01,
                             -1.3623e-01},
                            {6.2799e-02,
                             6.3892e-02,
                             6.4943e-02,
                             6.4881e-02,
                             7.2169e-02,
                             9.5669e-02,
                             1.3082e-01,
                             1.7694e-01,
                             2.2559e-01,
                             2.7655e-01,
                             3.2483e-01,
                             3.6888e-01,
                             4.0783e-01,
                             4.4330e-01,
                             4.7317e-01,
                             4.9631e-01,
                             5.1334e-01,
                             5.2318e-01,
                             5.2651e-01},
                            {-1.4949e+00,
                             -1.4952e+00,
                             -1.4925e+00,
                             -1.4889e+00,
                             -1.4867e+00,
                             -1.4897e+00,
                             -1.4948e+00,
                             -1.5025e+00,
                             -1.5084e+00,
                             -1.5142e+00,
                             -1.5176e+00,
                             -1.5191e+00,
                             -1.5187e+00,
                             -1.5180e+00,
                             -1.5169e+00,
                             -1.5153e+00,
                             -1.5138e+00,
                             -1.5125e+00,
                             -1.5120e+00},
                            {0.0000e+00,
                             0.0000e+00,
                             0.0000e+00,
                             0.0000e+00,
                             0.0000e+00,
                             0.0000e+00,
                             0.0000e+00,
                             0.0000e+00,
                             0.0000e+00,
                             0.0000e+00,
                             0.0000e+00,
                             0.0000e+00,
                             0.0000e+00,
                             0.0000e+00,
                             0.0000e+00,
                             0.0000e+00,
                             0.0000e+00,
                             0.0000e+00,
                             0.0000e+00}};

  double tof_pixel_spherical_absorption_correction(double muR,
                                                   double two_theta,
                                                   int two_theta_idx) {
    /*
     * Calculate the spherical absorption correction for a single pixel
     * defined by two_theta
     */

    double ln_t1 = 0;
    double ln_t2 = 0;
    for (std::size_t k = 0; k < 8; ++k) {
      ln_t1 = ln_t1 * muR + pc[k][two_theta_idx];
      ln_t2 = ln_t2 * muR + pc[k][two_theta_idx + 1];
    }
    const double t1 = exp(ln_t1);
    const double t2 = exp(ln_t2);
    const double sin_theta_1 = pow(sin(deg_as_rad(two_theta_idx * 5.0)), 2);
    const double sin_theta_2 = pow(sin(deg_as_rad((two_theta_idx + 1) * 5.0)), 2);
    const double l1 = (t1 - t2) / (sin_theta_1 - sin_theta_2);
    const double l0 = t1 - l1 * sin_theta_1;
    const double correction = 1 / (l0 + l1 * pow(sin(two_theta * .5), 2));
    return correction;
  }

  /*
   * Holds constants required for correcting time-of-flight data for absorption
   */
  struct TOFAbsorptionParams {
    double sample_radius;
    double sample_scattering_x_section;
    double sample_absorption_x_section;
    double sample_number_density;
    double incident_radius;
    double incident_scattering_x_section;
    double incident_absorption_x_section;
    double incident_number_density;
    double sample_linear_scattering_c;
    double incident_linear_scattering_c;
    double sample_linear_absorption_c;
    double incident_linear_absorption_c;

    TOFAbsorptionParams(double sample_radius,
                        double sample_scattering_x_section,
                        double sample_absorption_x_section,
                        double sample_number_density,
                        double incident_radius,
                        double incident_scattering_x_section,
                        double incident_absorption_x_section,
                        double incident_number_density)
        : sample_radius(sample_radius * .1),  // Given in mm but calculated in cm
          sample_scattering_x_section(sample_scattering_x_section),
          sample_absorption_x_section(sample_absorption_x_section),
          sample_number_density(sample_number_density),
          incident_radius(incident_radius * .1),  // Given in mm but calculated in cm
          incident_scattering_x_section(incident_scattering_x_section),
          incident_absorption_x_section(incident_absorption_x_section),
          incident_number_density(incident_number_density) {
      sample_linear_scattering_c = sample_number_density * sample_scattering_x_section;
      incident_linear_scattering_c =
        incident_number_density * incident_scattering_x_section;
      sample_linear_absorption_c = sample_number_density * sample_absorption_x_section;
      incident_linear_absorption_c =
        incident_number_density * incident_absorption_x_section;
    }
  };

  void extract_shoeboxes_to_reflection_table(
    dials::af::reflection_table &reflection_table,
    Experiment &experiment,
    ImageSequence &data) {
    Detector detector = *experiment.get_detector();
    Scan scan = *experiment.get_scan();

    scitbx::af::shared<double> img_tof = scan.get_property<double>("time_of_flight");

    int n_panels = detector.size();
    int num_images = data.size();
    vec2<std::size_t> image_size = detector[0].get_image_size();
    DIALS_ASSERT(num_images == img_tof.size());

    // Processor to get the image data into shoeboxes
    ShoeboxProcessor shoebox_processor(
      reflection_table, n_panels, 0, num_images, false);

    // Get shoeboxes for image data
    for (std::size_t img_num = 0; img_num < num_images; ++img_num) {
      dxtbx::format::Image<double> img = data.get_corrected_data(img_num);
      dxtbx::format::Image<bool> mask = data.get_mask(img_num);

      dials::af::shared<scitbx::af::versa<double, scitbx::af::c_grid<2>>> output_data(
        n_panels);
      dials::af::shared<scitbx::af::versa<bool, scitbx::af::c_grid<2>>> output_mask(
        n_panels);

      for (std::size_t i = 0; i < output_data.size(); ++i) {
        output_data[i] = img.tile(i).data();
        output_mask[i] = mask.tile(i).data();
      }
      shoebox_processor.next_data_only(
        dials::model::Image<double>(output_data.const_ref(), output_mask.const_ref()));
    }
  }

  void integrate_reflection_table(dials::af::reflection_table &reflection_table,
                                  Experiment &experiment,
                                  ImageSequence &data,
                                  bool apply_lorentz_correction) {
    /*
     * Updates reflection_table with intensities and variances with
     * optional Lorentz correction
     */

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

    dials::af::shared<bool> succeeded(reflection_table.size());
    dials::af::shared<double> intensities(reflection_table.size());
    dials::af::shared<double> variances(reflection_table.size());
    int bg_code = Valid | Background | BackgroundUsed;

    for (std::size_t i = 0; i < reflection_table.size(); ++i) {
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

        int frame_z = bbox[4] + z;
        double tof = img_tof[frame_z] * std::pow(10, -6);  // (s)

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
        if (!success) {
          break;
        }

        int frame_z = bbox[4] + z;
        double tof = img_tof[frame_z] * std::pow(10, -6);  // (s)

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
            double var_B =
              std::abs(raw_B) * (1.0 + double(n_signal) / double(n_background));

            // Lorentz correction
            scitbx::vec3<double> s1 = detector[panel].get_pixel_lab_coord(
              scitbx::vec2<double>(panel_x, panel_y));
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

            double I = I0 / bin_width_correction;
            double var_I = var_I0 / (bin_width_correction * bin_width_correction);

            if (apply_lorentz_correction) {
              I *= L;
              var_I *= (L * L);
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
      }
      succeeded[i] = success;
      intensities[i] = intensity;
      variances[i] = variance;
    }
    reflection_table["intensity.sum.value"] = intensities;
    reflection_table["intensity.sum.variance"] = variances;

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
    }
  }

  std::vector<double> savitzky_golay(const std::vector<double> &signal,
                                     int window_size,
                                     int poly_order) {
    /*
     * Method to smooth spectra to reduce noise
     * Based on https://en.wikipedia.org/wiki/Savitzky%E2%80%93Golay_filter
     */

    // window_size must be odd
    DIALS_ASSERT(window_size % 2 != 0);
    DIALS_ASSERT(window_size > 1);
    DIALS_ASSERT(poly_order < window_size);

    int half_window = (window_size - 1) / 2;
    int n = signal.size();
    std::vector<double> result(n, 0.0);

    // Precompute coefficients using least-squares polynomial fit
    std::vector<double> coeffs(window_size);
    {
      // Build Vandermonde matrix and compute pseudoinverse
      std::vector<std::vector<double>> A(window_size,
                                         std::vector<double>(poly_order + 1, 0.0));
      for (int i = -half_window; i <= half_window; ++i) {
        for (int j = 0; j <= poly_order; ++j) {
          A[i + half_window][j] = std::pow(i, j);
        }
      }

      // Compute pseudoinverse using normal equations (A^T A)^-1 A^T
      std::vector<std::vector<double>> ATA(poly_order + 1,
                                           std::vector<double>(poly_order + 1, 0.0));
      for (int i = 0; i <= poly_order; ++i) {
        for (int j = 0; j <= poly_order; ++j) {
          for (int k = 0; k < window_size; ++k) {
            ATA[i][j] += A[k][i] * A[k][j];
          }
        }
      }
      // Solve linear system ATA * c = b
      // where b = column corresponding to evaluating at center point
      std::vector<double> b(poly_order + 1, 0.0);
      b[0] = 1.0;  // evaluate polynomial at x=0

      // Gaussian elimination
      for (int i = 0; i <= poly_order; ++i) {
        double pivot = ATA[i][i];
        for (int j = i; j <= poly_order; ++j) {
          ATA[i][j] /= pivot;
        }
        b[i] /= pivot;
        for (int k = i + 1; k <= poly_order; ++k) {
          double factor = ATA[k][i];
          for (int j = i; j <= poly_order; ++j) {
            ATA[k][j] -= factor * ATA[i][j];
          }
          b[k] -= factor * b[i];
        }
      }
      for (int i = poly_order; i >= 0; --i) {
        for (int j = i + 1; j <= poly_order; ++j) {
          b[i] -= ATA[i][j] * b[j];
        }
      }
      // Compute convolution coeffs
      for (int k = -half_window; k <= half_window; ++k) {
        double sum = 0.0;
        for (int j = 0; j <= poly_order; ++j) {
          sum += b[j] * std::pow(k, j);
        }
        coeffs[k + half_window] = sum;
      }
    }

    // Convolve
    for (int i = 0; i < n; ++i) {
      double sm = 0.0;
      for (int j = -half_window; j <= half_window; ++j) {
        int idx = i + j;
        if (idx < 0) idx = 0;
        if (idx >= n) idx = n - 1;
        sm += coeffs[j + half_window] * signal[idx];
      }
      result[i] = sm;
    }

    return result;
  }

  void integrate_reflection_table(dials::af::reflection_table &reflection_table,
                                  Experiment &experiment,
                                  ImageSequence &data,
                                  ImageSequence &incident_data,
                                  ImageSequence &empty_data,
                                  double sample_proton_charge,
                                  double incident_proton_charge,
                                  double empty_proton_charge,
                                  bool apply_lorentz_correction) {
    /*
     * Updates reflection_table with intensities and variances corrected by
     * incident and empty runs, and an optional Lorentz correction
     */

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
      dxtbx::format::Image<double> img = incident_data.get_corrected_data(img_num);
      dxtbx::format::Image<bool> mask = incident_data.get_mask(img_num);

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
      dxtbx::format::Image<double> img = empty_data.get_corrected_data(img_num);
      dxtbx::format::Image<bool> mask = empty_data.get_mask(img_num);

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

    // Now correct and integrate each pixel for each shoebox
    dials::af::shared<Shoebox<>> shoeboxes = reflection_table["shoebox"];
    dials::af::const_ref<int6> bboxes = reflection_table["bbox"];

    dials::af::shared<bool> succeeded(reflection_table.size());
    dials::af::shared<double> intensities(reflection_table.size());
    dials::af::shared<double> variances(reflection_table.size());
    int bg_code = Valid | Background | BackgroundUsed;

    for (std::size_t i = 0; i < reflection_table.size(); ++i) {
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
      std::vector<double> incident_spectrum(shoebox.zsize(), 0.0);
      std::vector<double> empty_spectrum(shoebox.zsize(), 0.0);
      std::vector<std::size_t> n_contrib(shoebox.zsize(), 0);

      for (std::size_t z = 0; z < shoebox.zsize(); ++z) {
        if (!success) {
          break;
        }

        int frame_z = bbox[4] + z;
        double tof = img_tof[frame_z] * std::pow(10, -6);  // (s)

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
      std::vector<double> smoothed_empty = savitzky_golay(empty_spectrum, 7, 2);
      std::vector<double> smoothed_incident = savitzky_golay(incident_spectrum, 7, 2);

      // Second pass to integrate the shoeboxes
      for (std::size_t z = 0; z < shoebox.zsize(); ++z) {
        if (!success) {
          break;
        }

        int frame_z = bbox[4] + z;
        double tof = img_tof[frame_z] * std::pow(10, -6);  // (s)

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
            double S = raw_S / sample_proton_charge;
            double B = raw_B / sample_proton_charge;
            double J = raw_J / incident_proton_charge;
            double E = raw_E / empty_proton_charge;

            // Variances
            double var_S =
              std::abs(raw_S) / (sample_proton_charge * sample_proton_charge);
            double var_B =
              (std::abs(raw_B) / (sample_proton_charge * sample_proton_charge))
              * (1.0 + double(n_signal) / double(n_background));
            double var_J = (std::abs(raw_J) * n)
                           / (incident_proton_charge * incident_proton_charge * n * n);
            double var_E = (std::abs(raw_E) * n)
                           / (empty_proton_charge * empty_proton_charge * n * n);

            // Lorentz correction
            scitbx::vec3<double> s1 = detector[panel].get_pixel_lab_coord(
              scitbx::vec2<double>(panel_x, panel_y));
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

            // Apply Lorentz correction and normalise by incident spectrum
            double I, var_I;
            if (apply_lorentz_correction) {
              I = L * I0 / J;
              // Var( L*I0/J ) = (L/J)^2 Var(I0) + (L*I0/J^2)^2 Var(J)
              var_I = (L * L / (J * J)) * var_I0
                      + (L * L * I0 * I0 / (J * J * J * J)) * var_J;
            } else {
              I = I0 / J;
              // Var( I0/J ) = (1/J)^2 Var(I0) + (I0/J^2)^2 Var(J)
              var_I = (1 / (J * J)) * var_I0 + (I0 * I0 / (J * J * J * J)) * var_J;
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
      }
      succeeded[i] = success;
      intensities[i] = intensity;
      variances[i] = variance;
    }
    reflection_table["intensity.sum.value"] = intensities;
    reflection_table["intensity.sum.variance"] = variances;

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
    }
  }

  void integrate_reflection_table(dials::af::reflection_table &reflection_table,
                                  Experiment &experiment,
                                  ImageSequence &data,
                                  ImageSequence &incident_data,
                                  ImageSequence &empty_data,
                                  double sample_proton_charge,
                                  double incident_proton_charge,
                                  double empty_proton_charge,
                                  TOFAbsorptionParams &corrections_data,
                                  bool apply_lorentz_correction) {
    /*
     * Updates reflection_table with intensities and variances corrected by
     * incident and empty runs, spherical absorption, and an optional Lorentz
     * correction
     */

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
      dxtbx::format::Image<double> img = incident_data.get_corrected_data(img_num);
      dxtbx::format::Image<bool> mask = incident_data.get_mask(img_num);

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
      dxtbx::format::Image<double> img = empty_data.get_corrected_data(img_num);
      dxtbx::format::Image<bool> mask = empty_data.get_mask(img_num);

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

    // Now correct and integrate each pixel for each shoebox
    dials::af::shared<Shoebox<>> shoeboxes = reflection_table["shoebox"];
    dials::af::const_ref<int6> bboxes = reflection_table["bbox"];

    dials::af::shared<bool> succeeded(reflection_table.size());
    dials::af::shared<double> intensities(reflection_table.size());
    dials::af::shared<double> variances(reflection_table.size());
    int bg_code = Valid | Background | BackgroundUsed;

    for (std::size_t i = 0; i < reflection_table.size(); ++i) {
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
      std::vector<double> incident_spectrum(shoebox.zsize(), 0.0);
      std::vector<double> empty_spectrum(shoebox.zsize(), 0.0);
      std::vector<std::size_t> n_contrib(shoebox.zsize(), 0);

      for (std::size_t z = 0; z < shoebox.zsize(); ++z) {
        if (!success) {
          break;
        }

        int frame_z = bbox[4] + z;
        double tof = img_tof[frame_z] * std::pow(10, -6);  // (s)

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

            // Get required pixel data
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
      std::vector<double> smoothed_empty = savitzky_golay(empty_spectrum, 7, 2);
      std::vector<double> smoothed_incident = savitzky_golay(incident_spectrum, 7, 2);

      // Second pass to integrate the shoeboxes
      for (std::size_t z = 0; z < shoebox.zsize(); ++z) {
        if (!success) {
          break;
        }

        int frame_z = bbox[4] + z;
        double tof = img_tof[frame_z] * std::pow(10, -6);  // (s)

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

            // Pixel data will be divided by incident spectrum
            if (raw_J < 1e-7) {
              continue;
            }

            int mask = shoebox.mask(z, y, x);

            // Normalise by proton charge
            double S = raw_S / sample_proton_charge;
            double B = raw_B / sample_proton_charge;
            double J = raw_J / incident_proton_charge;
            double E = raw_E / empty_proton_charge;

            // Variances
            double var_S =
              std::abs(raw_S) / (sample_proton_charge * sample_proton_charge);
            double var_B =
              (std::abs(raw_B) / (sample_proton_charge * sample_proton_charge))
              * (1.0 + double(n_signal) / double(n_background));
            double var_J = (std::abs(raw_J) * n)
                           / (incident_proton_charge * incident_proton_charge * n * n);
            double var_E = (std::abs(raw_E) * n)
                           / (empty_proton_charge * empty_proton_charge * n * n);

            // Lorentz correction
            scitbx::vec3<double> s1 = detector[panel].get_pixel_lab_coord(
              scitbx::vec2<double>(panel_x, panel_y));
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
            double T = tof_pixel_spherical_absorption_correction(
              sample_muR, two_theta, two_theta_idx);

            // Pixel data will be divided by the absorption correction
            if (T < 1e-7) {
              continue;
            }
            double incident_muR =
              (corrections_data.incident_linear_scattering_c
               + (corrections_data.incident_linear_absorption_c / 1.8) * wl)
              * corrections_data.incident_radius;
            double J_T = tof_pixel_spherical_absorption_correction(
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

              var_I = (L * L / (J * J * T * T)) * var_I0
                      + (L * L * I0 * I0 / (J * J * J * J * T * T)) * var_J;
            } else {
              I = I0 / (J * T);
              var_I = (1 / (J * J * T * T)) * var_I0
                      + (I0 * I0 / (J * J * J * J * T * T)) * var_J;
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
      }
      succeeded[i] = success;
      intensities[i] = intensity;
      variances[i] = variance;
    }
    reflection_table["intensity.sum.value"] = intensities;
    reflection_table["intensity.sum.variance"] = variances;
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
    }
  }

  void integrate_reflection_table_profile1d(
    dials::af::reflection_table &reflection_table,
    Experiment &experiment,
    ImageSequence &data,
    bool apply_lorentz_correction,
    double A,
    double alpha,
    double beta,
    double sigma) {
    /*
     * Updates reflection_table with intensities and variances with
     * optional Lorentz correction
     */

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

    dials::af::shared<bool> succeeded(reflection_table.size());
    dials::af::shared<bool> succeeded_prf(reflection_table.size());
    dials::af::shared<double> intensities(reflection_table.size());
    dials::af::shared<double> variances(reflection_table.size());
    dials::af::shared<double> intensities_prf(reflection_table.size());
    dials::af::shared<double> variances_prf(reflection_table.size());
    int bg_code = Valid | Background | BackgroundUsed;

    for (std::size_t i = 0; i < reflection_table.size(); ++i) {
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

      std::vector<double> tof_z(shoebox.zsize(), 0.0);
      std::vector<double> projected_intensity(shoebox.zsize(), 0.0);
      std::vector<double> projected_variance(shoebox.zsize(), 0.0);
      // Second pass to perform actual integration
      for (std::size_t z = 0; z < shoebox.zsize(); ++z) {
        double intensity_z = 0;
        double variance_z = 0;

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

            int mask = shoebox.mask(z, y, x);

            // Variances
            double var_S = std::abs(raw_S);
            double var_B =
              std::abs(raw_B) * (1.0 + double(n_signal) / double(n_background));

            // Lorentz correction
            scitbx::vec3<double> s1 = detector[panel].get_pixel_lab_coord(
              scitbx::vec2<double>(panel_x, panel_y));
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

            double I = I0 / bin_width_correction;
            double var_I = var_I0 / (bin_width_correction * bin_width_correction);

            if (apply_lorentz_correction) {
              I *= L;
              var_I *= (L * L);
            }

            // Accumulate if pixel in foreground & valid
            if ((mask & Foreground) == Foreground && (mask & Valid) == Valid
                && (mask & Overlapped) == 0) {
              intensity_z += I;
              variance_z += var_I;
              intensity += I;
              variance += var_I;
            } else if ((mask & Foreground) == Foreground) {
              success = false;
              break;
            }
          }
          projected_intensity[z] = intensity_z;
          projected_variance[z] = variance_z;
        }
      }
      succeeded[i] = success;
      intensities[i] = intensity;
      variances[i] = variance;

      if (success) {
        auto max_it =
          std::max_element(projected_intensity.begin(), projected_intensity.end());
        size_t max_index = std::distance(projected_intensity.begin(), max_it);
        double T_ph = tof_z[max_index];
        TOFProfile1D profile(tof_z, projected_intensity, A, alpha, beta, sigma, T_ph);

        bool profile_success = profile.fit();
        if (profile_success) {
          double I_prf = profile.calc_intensity();
          double var_prf = profile.calc_variance(projected_variance);
          if (var_prf > 1e-7) {
            double i_sigma_sum = intensity / std::sqrt(variance);
            double i_sigma_prf = I_prf / std::sqrt(var_prf);
            if (i_sigma_prf > i_sigma_sum) {
              intensities_prf[i] = I_prf;
              variances_prf[i] = var_prf;
            } else {
              profile_success = false;
            }
          } else {
            profile_success = false;
          }
        }
        succeeded_prf[i] = profile_success;
      } else {
        succeeded_prf[i] = false;
      }
    }

    reflection_table["intensity.sum.value"] = intensities;
    reflection_table["intensity.sum.variance"] = variances;
    reflection_table["intensity.prf.value"] = intensities_prf;
    reflection_table["intensity.prf.variance"] = variances_prf;

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
      if (succeeded_prf[i]) {
        flags[i] &= ~dials::af::FailedDuringProfileFitting;
        flags[i] |= dials::af::IntegratedPrf;
      } else {
        flags[i] &= ~dials::af::IntegratedPrf;
        flags[i] |= dials::af::FailedDuringProfileFitting;
      }
    }
  }

  void integrate_reflection_table_profile1d(
    dials::af::reflection_table &reflection_table,
    Experiment &experiment,
    ImageSequence &data,
    ImageSequence &incident_data,
    ImageSequence &empty_data,
    double sample_proton_charge,
    double incident_proton_charge,
    double empty_proton_charge,
    bool apply_lorentz_correction,
    double A,
    double alpha,
    double beta,
    double sigma) {
    /*
     * Updates reflection_table with intensities and variances corrected by
     * incident and empty runs, and an optional Lorentz correction
     */

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
      dxtbx::format::Image<double> img = incident_data.get_corrected_data(img_num);
      dxtbx::format::Image<bool> mask = incident_data.get_mask(img_num);

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
      dxtbx::format::Image<double> img = empty_data.get_corrected_data(img_num);
      dxtbx::format::Image<bool> mask = empty_data.get_mask(img_num);

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

    dials::af::shared<bool> succeeded(reflection_table.size());
    dials::af::shared<bool> succeeded_prf(reflection_table.size());
    dials::af::shared<double> intensities(reflection_table.size());
    dials::af::shared<double> variances(reflection_table.size());
    dials::af::shared<double> intensities_prf(reflection_table.size());
    dials::af::shared<double> variances_prf(reflection_table.size());
    int bg_code = Valid | Background | BackgroundUsed;

    for (std::size_t i = 0; i < reflection_table.size(); ++i) {
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
      std::vector<double> incident_spectrum(shoebox.zsize(), 0.0);
      std::vector<double> empty_spectrum(shoebox.zsize(), 0.0);
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
      std::vector<double> smoothed_empty = savitzky_golay(empty_spectrum, 7, 2);
      std::vector<double> smoothed_incident = savitzky_golay(incident_spectrum, 7, 2);

      std::vector<double> tof_z(shoebox.zsize(), 0.0);
      std::vector<double> projected_intensity(shoebox.zsize(), 0.0);
      std::vector<double> projected_variance(shoebox.zsize(), 0.0);
      // Second pass to perform actual integration
      for (std::size_t z = 0; z < shoebox.zsize(); ++z) {
        double intensity_z = 0;
        double variance_z = 0;

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

            // Normalise by proton charge
            double S = raw_S / sample_proton_charge;
            double B = raw_B / sample_proton_charge;
            double J = raw_J / incident_proton_charge;
            double E = raw_E / empty_proton_charge;

            // Variances
            double var_S =
              std::abs(raw_S) / (sample_proton_charge * sample_proton_charge);
            double var_B =
              (std::abs(raw_B) / (sample_proton_charge * sample_proton_charge))
              * (1.0 + double(n_signal) / double(n_background));
            double var_J = (std::abs(raw_J) * n)
                           / (incident_proton_charge * incident_proton_charge * n * n);
            double var_E = (std::abs(raw_E) * n)
                           / (empty_proton_charge * empty_proton_charge * n * n);

            // Lorentz correction
            scitbx::vec3<double> s1 = detector[panel].get_pixel_lab_coord(
              scitbx::vec2<double>(panel_x, panel_y));
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

            // Apply Lorentz correction and normalise by incident spectrum
            double I, var_I;
            if (apply_lorentz_correction) {
              I = L * I0 / J;
              // Var( L*I0/J ) = (L/J)^2 Var(I0) + (L*I0/J^2)^2 Var(J)
              var_I = (L * L / (J * J)) * var_I0
                      + (L * L * I0 * I0 / (J * J * J * J)) * var_J;
            } else {
              I = I0 / J;
              // Var( I0/J ) = (1/J)^2 Var(I0) + (I0/J^2)^2 Var(J)
              var_I = (1 / (J * J)) * var_I0 + (I0 * I0 / (J * J * J * J)) * var_J;
            }

            // Accumulate if pixel in foreground & valid
            if ((mask & Foreground) == Foreground && (mask & Valid) == Valid
                && (mask & Overlapped) == 0) {
              intensity_z += I;
              variance_z += var_I;
              intensity += I;
              variance += var_I;
            } else if ((mask & Foreground) == Foreground) {
              success = false;
              break;
            }
          }
          projected_intensity[z] = intensity_z;
          projected_variance[z] = variance_z;
        }
      }
      succeeded[i] = success;
      intensities[i] = intensity;
      variances[i] = variance;

      if (success) {
        auto max_it =
          std::max_element(projected_intensity.begin(), projected_intensity.end());
        size_t max_index = std::distance(projected_intensity.begin(), max_it);
        double T_ph = tof_z[max_index];
        TOFProfile1D profile(tof_z, projected_intensity, A, alpha, beta, sigma, T_ph);

        bool profile_success = profile.fit();
        if (profile_success) {
          double I_prf = profile.calc_intensity();
          double var_prf = profile.calc_variance(projected_variance);
          if (var_prf > 1e-7) {
            double i_sigma_sum = intensity / std::sqrt(variance);
            double i_sigma_prf = I_prf / std::sqrt(var_prf);
            if (i_sigma_prf > i_sigma_sum) {
              intensities_prf[i] = I_prf;
              variances_prf[i] = var_prf;
            } else {
              profile_success = false;
            }
          } else {
            profile_success = false;
          }
        }
        succeeded_prf[i] = profile_success;
      } else {
        succeeded_prf[i] = false;
      }
    }

    reflection_table["intensity.sum.value"] = intensities;
    reflection_table["intensity.sum.variance"] = variances;
    reflection_table["intensity.prf.value"] = intensities_prf;
    reflection_table["intensity.prf.variance"] = variances_prf;

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
      if (succeeded_prf[i]) {
        flags[i] &= ~dials::af::FailedDuringProfileFitting;
        flags[i] |= dials::af::IntegratedPrf;
      } else {
        flags[i] &= ~dials::af::IntegratedPrf;
        flags[i] |= dials::af::FailedDuringProfileFitting;
      }
    }
  }

  void integrate_shoebox_profile1d(Shoebox<> &shoebox,
                                   Experiment &experiment,
                                   ImageSequence &data,
                                   bool apply_lorentz_correction,
                                   double A,
                                   double alpha,
                                   double beta,
                                   double sigma,
                                   dials::af::shared<double> projected_intensity_out,
                                   dials::af::shared<double> line_profile_out,
                                   double &prf_intensity_out,
                                   double &prf_variance_out,
                                   double &sum_intensity_out,
                                   double &sum_variance_out,
                                   bool &profile_success) {
    /*
     * Calculates projected_intensity, line_profile, prf_intensity,
     * prf_variance, sum_intensity, sum_variance
     */

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

    int bg_code = Valid | Background | BackgroundUsed;

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

    std::vector<double> tof_z(shoebox.zsize(), 0.0);
    std::vector<double> projected_intensity(shoebox.zsize(), 0.0);
    std::vector<double> projected_variance(shoebox.zsize(), 0.0);
    // Second pass to perform actual integration
    for (std::size_t z = 0; z < shoebox.zsize(); ++z) {
      double intensity_z = 0;
      double variance_z = 0;

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

          int mask = shoebox.mask(z, y, x);

          // Variances
          double var_S = std::abs(raw_S);
          double var_B =
            std::abs(raw_B) * (1.0 + double(n_signal) / double(n_background));

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

          double I = I0 / bin_width_correction;
          double var_I = var_I0 / (bin_width_correction * bin_width_correction);

          if (apply_lorentz_correction) {
            I *= L;
            var_I *= (L * L);
          }

          // Accumulate if pixel in foreground & valid
          if ((mask & Foreground) == Foreground && (mask & Valid) == Valid
              && (mask & Overlapped) == 0) {
            intensity_z += I;
            variance_z += var_I;
            intensity += I;
            variance += var_I;
          } else if ((mask & Foreground) == Foreground) {
            success = false;
            break;
          }
        }
        projected_intensity[z] = intensity_z;
        projected_variance[z] = variance_z;
      }
    }

    if (success) {
      auto max_it =
        std::max_element(projected_intensity.begin(), projected_intensity.end());
      size_t max_index = std::distance(projected_intensity.begin(), max_it);
      double T_ph = tof_z[max_index];
      TOFProfile1D profile(tof_z, projected_intensity, A, alpha, beta, sigma, T_ph);

      profile_success = profile.fit();
      if (profile_success) {
        double I_prf = profile.calc_intensity();
        double var_prf = profile.calc_variance(projected_variance);
        if (var_prf > 1e-7) {
          double i_sigma_sum = intensity / std::sqrt(variance);
          double i_sigma_prf = I_prf / std::sqrt(var_prf);
          if (i_sigma_prf > i_sigma_sum) {
            prf_intensity_out = I_prf;
            prf_variance_out = var_prf;
            sum_intensity_out = intensity;
            sum_variance_out = variance;
            std::vector<double> line_profile = profile.result();
            DIALS_ASSERT(line_profile.size() == line_profile_out.size());
            DIALS_ASSERT(projected_intensity.size() == projected_intensity_out.size());
            for (std::size_t i = 0; i < line_profile.size(); ++i) {
              line_profile_out[i] = line_profile[i];
              projected_intensity_out[i] = projected_intensity[i];
            }
          } else {
            profile_success = false;
          }
        } else {
          profile_success = false;
        }
      }
    } else {
      profile_success = false;
    }
  }

}}  // namespace dials::algorithms

#endif /* DIALS_ALGORITHMS_INTEGRATION_TOF_INTEGRATION_H */