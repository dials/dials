#ifndef DIALS_ALGORITHMS_INTEGRATION_TOF_MASK_CALCULATOR_H
#define DIALS_ALGORITHMS_INTEGRATION_TOF_MASK_CALCULATOR_H

#include <dxtbx/imageset.h>
#include <dxtbx/format/image.h>
#include <dxtbx/array_family/flex_table.h>
#include <dials/array_family/reflection_table.h>
#include <dials/model/data/shoebox.h>
#include <dxtbx/model/detector.h>
#include <dxtbx/model/beam.h>
#include <dxtbx/model/scan.h>
#include <dxtbx/model/experiment.h>
#include <dxtbx/model/goniometer.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include <scitbx/constants.h>
#include <eigen3/Eigen/Dense>
#include <vector>

namespace dials { namespace algorithms {

  using dials::model::Shoebox;
  using dxtbx::ImageSequence;
  using dxtbx::af::flex_table;
  using dxtbx::model::Detector;
  using dxtbx::model::Experiment;
  using dxtbx::model::Goniometer;
  using dxtbx::model::PolychromaticBeam;
  using dxtbx::model::Scan;
  using dxtbx::model::scan_property_types;
  using scitbx::mat3;
  using scitbx::vec2;
  using scitbx::vec3;
  using scitbx::af::int6;
  using scitbx::constants::m_n;
  using scitbx::constants::pi;
  using scitbx::constants::Planck;

  void compute_weighted_ellipsoid(std::vector<Eigen::Vector3d>& points,
                                  std::vector<double>& values,
                                  Eigen::Vector3d& mean,
                                  Eigen::Matrix3d& eigenvectors,
                                  Eigen::Vector3d& axes_lengths,
                                  bool calculate_mean) {
    if (points.size() != values.size() || points.empty()) {
      throw std::invalid_argument(
        "Points and values must have the same size and cannot be empty.");
    }

    // Normalize values to compute weights
    double totalWeight = std::accumulate(values.begin(), values.end(), 0.0);
    if (totalWeight == 0) {
      throw std::invalid_argument("Values must not all be zero.");
    }

    std::vector<double> weights(values.size());
    for (size_t i = 0; i < values.size(); ++i) {
      weights[i] = values[i] / totalWeight;
    }

    if (calculate_mean) {
      // Compute the weighted mean
      Eigen::Vector3d weightedSum = Eigen::Vector3d::Zero();
      for (size_t i = 0; i < points.size(); ++i) {
        weightedSum += weights[i] * points[i];
      }

      mean = weightedSum;
    }

    // Compute the weighted covariance matrix
    Eigen::Matrix3d covMatrix = Eigen::Matrix3d::Zero();
    for (size_t i = 0; i < points.size(); ++i) {
      Eigen::Vector3d centered = points[i] - mean;
      covMatrix += weights[i] * (centered * centered.transpose());
    }

    // Perform eigen decomposition on the covariance matrix
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(covMatrix);
    if (eigensolver.info() != Eigen::Success) {
      throw std::runtime_error("Eigen decomposition failed.");
    }

    // Eigenvalues and eigenvectors
    Eigen::VectorXd eigenvalues = eigensolver.eigenvalues();
    eigenvectors = eigensolver.eigenvectors();
    axes_lengths =
      eigenvalues.cwiseSqrt();  // Semi-axis lengths are sqrt of eigenvalues
  }

  bool point_inside_ellipsoid(const Eigen::Vector3d& point,
                              const Eigen::Vector3d& mean,
                              const Eigen::Matrix3d& eigenvectors,
                              const Eigen::Vector3d& axes_lengths) {
    Eigen::Vector3d centered_point = point - mean;
    Eigen::Vector3d transformed_point = eigenvectors.transpose() * centered_point;
    Eigen::Vector3d normalized_point = transformed_point.cwiseQuotient(axes_lengths);
    double distance_squared = normalized_point.squaredNorm();
    return distance_squared <= 1.0;
  }

  void get_shoebox_rlps(const Shoebox<>& shoebox,
                        std::vector<Eigen::Vector3d>& shoebox_rlps,
                        std::vector<double>& shoebox_values,
                        Detector& detector,
                        int& panel,
                        int6& bbox,
                        scitbx::af::shared<double>& img_tof,
                        vec3<double>& unit_s0,
                        double& sample_to_source_distance,
                        mat3<double>& setting_rotation) {
    /*
     * Fills shoebox_rlps with beam vectors for each pixel in shoebox
     * Fills shoebox_values with the data of shoebox.data
     */

    for (std::size_t z = 0; z < shoebox.zsize(); ++z) {
      int frame_z = bbox[4] + z;
      double tof = img_tof[frame_z] * std::pow(10, -6);  // (s)

      for (std::size_t y = 0; y < shoebox.ysize(); ++y) {
        int panel_y = bbox[2] + y;

        for (std::size_t x = 0; x < shoebox.xsize(); ++x) {
          int panel_x = bbox[0] + x;
          vec3<double> s1 =
            detector[panel].get_pixel_lab_coord(vec2<double>(panel_x, panel_y));
          double pixel_distance = s1.length() + sample_to_source_distance;
          pixel_distance *= std::pow(10, -3);  // (m)
          double wl =
            ((Planck * tof) / (m_n * (pixel_distance))) * std::pow(10, 10);  // (A)
          vec3<double> s0 = unit_s0 / wl;
          s1 = s1 / s1.length() * (1 / wl);
          vec3<double> S = s1 - s0;
          S = setting_rotation.inverse() * S;
          shoebox_rlps.emplace_back(Eigen::Vector3d(S[0], S[1], S[2]));
          shoebox_values.push_back(shoebox.data(z, y, x));
        }
      }
    }
  }

  std::vector<size_t> get_shoebox_pixel_neighbors(size_t index,
                                                  std::size_t zsize,
                                                  std::size_t ysize,
                                                  std::size_t xsize) {
    std::vector<size_t> neighbors;

    std::size_t z = index / (ysize * xsize);
    std::size_t y = (index % (ysize * xsize)) / xsize;
    std::size_t x = index % xsize;

    for (int dz = -1; dz <= 1; ++dz) {
      for (int dy = -1; dy <= 1; ++dy) {
        for (int dx = -1; dx <= 1; ++dx) {
          if (dz == 0 && dy == 0 && dx == 0) continue;  // Skip the center pixel
          int nz = z + dz;
          int ny = y + dy;
          int nx = x + dx;
          if (nz >= 0 && nz < zsize && ny >= 0 && ny < ysize && nx >= 0 && nx < xsize) {
            neighbors.push_back(nz * (ysize * xsize) + ny * xsize + nx);
          }
        }
      }
    }
    return neighbors;
  }

  double calculate_skewness(af::ref<float, af::c_grid<3>> data,
                            std::set<size_t> selected_pixels) {
    float avg_intensity = 0;
    int count = 0;
    for (std::size_t i = 0; i < data.size(); ++i) {
      if (selected_pixels.find(i) == selected_pixels.end()) {
        avg_intensity += data[i];
        count++;
      }
    }
    avg_intensity /= count;

    float skewness = 0;
    for (std::size_t i = 0; i < data.size(); ++i) {
      if (selected_pixels.find(i) == selected_pixels.end()) {
        skewness += std::pow(data[i] - avg_intensity, 3);
      }
    }
    skewness /= count;
    return skewness;
  }

  af::ref<int, af::c_grid<3>> fill_holes_in_mask(af::ref<int, af::c_grid<3>>& mask) {
    /*
     * Flood fill to ensure no holes in mask
     */

    std::size_t zsize = mask.accessor()[0];
    std::size_t ysize = mask.accessor()[1];
    std::size_t xsize = mask.accessor()[2];

    auto is_valid = [&](int z, int y, int x) {
      return z >= 0 && z < zsize && y >= 0 && y < ysize && x >= 0 && x < xsize;
    };

    auto is_background = [&](int z, int y, int x) {
      return (mask(z, y, x) & Background) == Background;
    };

    auto is_foreground = [&](int z, int y, int x) {
      return (mask(z, y, x) & Foreground) == Foreground;
    };

    af::versa<bool, af::c_grid<3>> visited(mask.accessor(), false);

    auto flood_fill = [&](
                        int start_z, int start_y, int start_x, bool& touches_boundary) {
      std::vector<std::tuple<int, int, int>> stack;   // Stack for flood-fill
      std::vector<std::tuple<int, int, int>> region;  // To store pixels in the region
      stack.emplace_back(start_z, start_y, start_x);

      // Local visited map for the current flood-fill
      af::versa<bool, af::c_grid<3>> fill_visited(visited.accessor(), false);
      std::copy(visited.begin(), visited.end(), fill_visited.begin());

      while (!stack.empty()) {
        auto tuple = stack.back();
        stack.pop_back();

        int z = std::get<0>(tuple);
        int y = std::get<1>(tuple);
        int x = std::get<2>(tuple);

        if (z == 0 || z == zsize - 1 || y == 0 || y == ysize - 1 || x == 0
            || x == xsize - 1) {
          touches_boundary = true;
          break;
        }

        if (!is_valid(z, y, x) || fill_visited(z, y, x) || is_foreground(z, y, x)) {
          continue;
        }

        fill_visited(z, y, x) = true;  // Temporarily mark as visited in this flood-fill
        region.emplace_back(z, y, x);

        // Add neighbors
        for (int dz = -1; dz <= 1; ++dz) {
          for (int dy = -1; dy <= 1; ++dy) {
            for (int dx = -1; dx <= 1; ++dx) {
              if ((std::abs(dz) + std::abs(dy) + std::abs(dx)) == 1) {
                int nz = z + dz, ny = y + dy, nx = x + dx;
                if (is_valid(nz, ny, nx) && !fill_visited(nz, ny, nx)) {
                  stack.emplace_back(nz, ny, nx);
                }
              }
            }
          }
        }
      }

      // If the region does not touch the boundary, update the global visited map
      if (!touches_boundary) {
        for (const auto& pixel : region) {
          int rz = std::get<0>(pixel);
          int ry = std::get<1>(pixel);
          int rx = std::get<2>(pixel);
          visited(rz, ry, rx) = true;
        }
      }

      return region;
    };

    // Loop over all pixels in the mask
    int count = 0;
    for (std::size_t z = 0; z < zsize; ++z) {
      for (std::size_t y = 0; y < ysize; ++y) {
        for (std::size_t x = 0; x < xsize; ++x) {
          // Look for unvisited Background pixels
          if (!visited(z, y, x) && is_background(z, y, x)) {
            bool touches_boundary = false;
            auto region = flood_fill(z, y, x, touches_boundary);

            // Only fill regions that do not touch the boundary
            if (!touches_boundary) {
              for (const auto& pixel : region) {
                int rz = std::get<0>(pixel);
                int ry = std::get<1>(pixel);
                int rx = std::get<2>(pixel);

                mask(rz, ry, rx) &= ~Background;
                mask(rz, ry, rx) |= Foreground;
                count++;
              }
            }
          }
        }
      }
    }

    return mask;
  }

  void tof_calculate_seed_skewness_shoebox_mask(af::reflection_table& reflection_table,
                                                Experiment& experiment,
                                                float d_skewness_threshold,
                                                int min_iterations) {
    /**
     * Implementation based on
     * J. Peters, The 'seed-skewness' integration method generalized for
     * three-dimensional Bragg peaks, (2003), 36, 1475-1479
     * https://doi.org/10.1107/S0021889803021939
     *
     * The shoebox masks in reflection_table are updated based on this algorithm
     *
     * @param d_skewness_threshold When change in skewness falls below this value
     * the algorithm stops
     * @param min_iterations: Minimum number of iterations before checking against
     * d_skewness_threshold for convergence
     */

    af::shared<Shoebox<>> shoeboxes = reflection_table["shoebox"];

    for (std::size_t i = 0; i < reflection_table.size(); ++i) {
      // Get shoebox data
      Shoebox<> shoebox = shoeboxes[i];
      af::ref<int, af::c_grid<3>> mask = shoebox.mask.ref();
      af::ref<float, af::c_grid<3>> data = shoebox.data.ref();
      std::size_t zsize = data.accessor()[0];
      std::size_t ysize = data.accessor()[1];
      std::size_t xsize = data.accessor()[2];

      // Keep track of which pixels are selected
      std::set<std::size_t> selected_pixels;

      // Initial seed
      auto max_it = std::max_element(data.begin(), data.end());
      DIALS_ASSERT(max_it != data.end());  // data empty
      float max_val = *max_it;
      std::size_t max_idx = std::distance(data.begin(), max_it);
      selected_pixels.insert(max_idx);

      float skewness = calculate_skewness(data, selected_pixels);
      float d_skewness = 0.0;
      int num_iterations = 0;

      while ((d_skewness < 0 && std::abs(d_skewness) > d_skewness_threshold)
             || (num_iterations < min_iterations)) {
        std::set<std::size_t> neighbors;

        // Collect neighbors of all selected pixels
        for (std::size_t idx : selected_pixels) {
          auto neighbor_indices = get_shoebox_pixel_neighbors(idx, zsize, ysize, xsize);
          for (std::size_t n_idx : neighbor_indices) {
            neighbors.insert(n_idx);
          }
        }

        // Find the neighbor with the maximum value not in selected pixels
        auto max_it = std::max_element(
          neighbors.begin(), neighbors.end(), [&](std::size_t a, std::size_t b) {
            bool a_selected = selected_pixels.find(a) != selected_pixels.end();
            bool b_selected = selected_pixels.find(b) != selected_pixels.end();

            if (a_selected && b_selected) return false;
            if (a_selected) return true;
            if (b_selected) return false;

            return data[a] < data[b];
          });

        if (max_it != neighbors.end()
            && selected_pixels.find(*max_it) == selected_pixels.end()) {
          std::size_t max_neighbor_index = *max_it;
          float max_neighbor_value = data[*max_it];
          selected_pixels.insert(max_neighbor_index);
        } else {
          break;
        }

        float new_skewness = calculate_skewness(data, selected_pixels);
        d_skewness = new_skewness - skewness;
        skewness = new_skewness;
        num_iterations++;
      }

      // Now update the shoebox mask with selected pixels as foreground
      for (std::size_t z = 0; z < zsize; ++z) {
        for (std::size_t y = 0; y < ysize; ++y) {
          for (std::size_t x = 0; x < xsize; ++x) {
            std::size_t idx = z * (ysize * xsize) + y * xsize + x;
            int mask_value = selected_pixels.find(idx) != selected_pixels.end()
                               ? Foreground
                               : Background;
            mask(z, y, x) &= ~(Foreground | Background);
            mask(z, y, x) |= mask_value;
          }
        }
      }

      mask = fill_holes_in_mask(mask);
    }
  }

  void tof_calculate_ellipse_shoebox_mask(af::reflection_table& reflection_table,
                                          Experiment& experiment) {
    /**
     * Updates the masks of shoeboxes in reflection_table based on weighted
     * ellipses in reciprocal space
     */

    af::shared<Shoebox<>> shoeboxes = reflection_table["shoebox"];
    Scan scan = *experiment.get_scan();
    Detector detector = *experiment.get_detector();

    // Set setting_rotation to identity if no goniometer present
    mat3<double> setting_rotation;
    std::shared_ptr<Goniometer> goniometer_ptr = experiment.get_goniometer();
    if (goniometer_ptr == nullptr) {
      setting_rotation = mat3<double>(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
    } else {
      setting_rotation = goniometer_ptr->get_setting_rotation();
    }

    // Ensure beam is polychromatic
    std::shared_ptr<dxtbx::model::BeamBase> beam_ptr = experiment.get_beam();
    std::shared_ptr<PolychromaticBeam> beam =
      std::dynamic_pointer_cast<PolychromaticBeam>(beam_ptr);
    DIALS_ASSERT(beam != nullptr);
    vec3<double> unit_s0 = beam->get_unit_s0();
    double sample_to_source_distance = beam->get_sample_to_source_distance();

    scitbx::af::shared<double> img_tof = scan.get_property<double>("time_of_flight");
    af::const_ref<int6> bboxes = reflection_table["bbox"];
    scitbx::af::shared<vec3<double>> rlps = reflection_table["rlp"];

    for (std::size_t i = 0; i < reflection_table.size(); ++i) {
      /*
       * For each reflection get rlps for each shoebox pixel in detector space
       * Then compute an ellipsoid based on the intensity of these points and
       * check which points sit inside it
       */

      Shoebox<> shoebox = shoeboxes[i];
      af::ref<int, af::c_grid<3>> mask = shoebox.mask.ref();
      int panel = shoebox.panel;
      int6 bbox = bboxes[i];
      vec3<double> rlp = rlps[i];
      std::vector<Eigen::Vector3d> shoebox_rlps;
      std::vector<double> shoebox_values;
      get_shoebox_rlps(shoebox,
                       shoebox_rlps,
                       shoebox_values,
                       detector,
                       panel,
                       bbox,
                       img_tof,
                       unit_s0,
                       sample_to_source_distance,
                       setting_rotation);

      // Centre the ellipse around the peak value
      auto peak_val = std::max_element(shoebox_values.begin(), shoebox_values.end());
      size_t peak_idx = std::distance(shoebox_values.begin(), peak_val);
      Eigen::Vector3d mean = shoebox_rlps[peak_idx];
      Eigen::Matrix3d eigenvectors;
      Eigen::Vector3d axes_lengths;
      compute_weighted_ellipsoid(
        shoebox_rlps, shoebox_values, mean, eigenvectors, axes_lengths, false);
      int count = 0;
      for (std::size_t z = 0; z < shoebox.zsize(); ++z) {
        for (std::size_t y = 0; y < shoebox.ysize(); ++y) {
          for (std::size_t x = 0; x < shoebox.xsize(); ++x) {
            int mask_value = point_inside_ellipsoid(
                               shoebox_rlps[count], mean, eigenvectors, axes_lengths)
                               ? Foreground
                               : Background;
            mask(z, y, x) &= ~(Foreground | Background);
            mask(z, y, x) |= mask_value;
            count++;
          }
        }
      }
    }
  }

}}  // namespace dials::algorithms

#endif /* DIALS_ALGORITHMS_INTEGRATION_TOF_MASK_CALCULATOR_H */