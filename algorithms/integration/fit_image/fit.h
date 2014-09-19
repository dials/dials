/*
 * fit.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */

#ifndef DIALS_ALGORITHMS_INTEGRATION_FIT_IMAGE_FIT_H
#define DIALS_ALGORITHMS_INTEGRATION_FIT_IMAGE_FIT_H

#include <scitbx/array_family/tiny_types.h>
#include <dxtbx/model/beam.h>
#include <dxtbx/model/detector.h>
#include <dxtbx/model/goniometer.h>
#include <dxtbx/model/scan.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/array_family/reflection_table.h>
#include <dials/algorithms/integration/profile/fitting.h>
#include <dials/algorithms/reflection_basis/coordinate_system.h>
#include <dials/algorithms/polygon/clip/clip.h>
#include <dials/algorithms/polygon/spatial_interpolation.h>
#include <dials/algorithms/polygon/area.h>


namespace dials { namespace algorithms {

  using scitbx::af::int3;
  using dxtbx::model::Beam;
  using dxtbx::model::Detector;
  using dxtbx::model::Goniometer;
  using dxtbx::model::Scan;
  using dxtbx::model::Panel;
  using dials::af::ReferenceSpot;
  using dials::af::DontIntegrate;
  using dials::af::IntegratedPrf;
  using dials::model::Shoebox;
  using dials::model::Valid;
  using dials::model::Foreground;
  using dials::algorithms::reflection_basis::CoordinateSystem;
  using dials::algorithms::polygon::simple_area;
  using dials::algorithms::polygon::clip::vert4;
  using dials::algorithms::polygon::clip::vert8;
  using dials::algorithms::polygon::clip::quad_with_convex_quad;
  using dials::algorithms::polygon::spatial_interpolation::reverse_quad_inplace_if_backward;

  template <typename T>
  T min4(T a, T b, T c, T d) {
    return std::min(std::min(a,b), std::min(c,d));
  }

  template <typename T>
  T max4(T a, T b, T c, T d) {
    return std::max(std::max(a,b), std::max(c,d));
  }

  /**
   * A class to specify the experiments input
   */
  class Spec {
  public:

    Spec(const Beam &beam,
         const Detector &detector,
         const Goniometer &goniometer,
         const Scan &scan,
         double delta_b,
         double delta_m)
      : beam_(beam),
        detector_(detector),
        goniometer_(goniometer),
        scan_(scan),
        delta_b_(delta_b),
        delta_m_(delta_m) {
      DIALS_ASSERT(delta_b_ > 0);
      DIALS_ASSERT(delta_m_ > 0);
    }

    const Beam& beam() const {
      return beam_;
    }

    const Detector& detector() const {
      return detector_;
    }

    const Goniometer& goniometer() const {
      return goniometer_;
    }

    const Scan& scan() const {
      return scan_;
    }

    const double delta_b() const {
      return delta_b_;
    }

    const double delta_m() const {
      return delta_m_;
    }

  private:

    Beam beam_;
    Detector detector_;
    Goniometer goniometer_;
    Scan scan_;
    double delta_b_;
    double delta_m_;
  };


  /**
   * A class to compute a single reference profile.
   */
  template <typename FloatType>
  class SingleProfileLearner {
  public:

    typedef FloatType float_type;
    typedef af::versa< FloatType, af::c_grid<3> > profile_type;
    typedef af::versa< int, af::c_grid<3> > mask_type;

    /**
     * Initialise the reference learner.
     * @param spec The experiment spec
     * @param grid_size The size of the reference grid
     */
    SingleProfileLearner(const Spec spec,
                         int3 grid_size)
        : detector_(spec.detector()),
          scan_(spec.scan()),
          m2_(spec.goniometer().get_rotation_axis()),
          s0_(spec.beam().get_s0()),
          delta_b_(spec.delta_b()),
          delta_m_(spec.delta_m()),
          data_(af::c_grid<3>(
                2 * grid_size[0] + 1,
                2 * grid_size[1] + 1,
                2 * grid_size[2] + 1), 0.0),
          finalized_(false),
          count_(0) {
      DIALS_ASSERT(grid_size[0] > 0);
      DIALS_ASSERT(grid_size[1] > 0);
      DIALS_ASSERT(grid_size[2] > 0);
    }

    /**
     * Add a reflection to the reference profile.
     * @param sbox The shoebox
     * @param s1 The diffracted beam vector.
     * @param phi The rotation angle.
     */
    bool add(const Shoebox<> &sbox, vec3<double> s1, double phi) {

      // FIXME Offset along e2 axis

      // Check the input
      DIALS_ASSERT(sbox.is_consistent());

      // Get the bbox size
      std::size_t xs = sbox.xsize();
      std::size_t ys = sbox.ysize();
      std::size_t zs = sbox.zsize();

      // Compute the grid step and offset
      double xoff = -delta_b_;
      double yoff = -delta_b_;
      double zoff = -delta_m_;
      double xstep = (2.0 * delta_b_) / data_.accessor()[2];
      double ystep = (2.0 * delta_b_) / data_.accessor()[1];
      double zstep = (2.0 * delta_m_) / data_.accessor()[0];

      // Create the coordinate system
      CoordinateSystem cs(m2_, s0_, s1, phi);

      // Get the panel
      const Panel &panel = detector_[sbox.panel];

      // Compute the detector coordinates of each point on the grid
      af::versa< vec2<double>, af::c_grid<2> > xy(af::c_grid<2>(
            data_.accessor()[1] + 1,
            data_.accessor()[2] + 1));
      for (std::size_t j = 0; j <= data_.accessor()[1]; ++j) {
        for (std::size_t i = 0; i <= data_.accessor()[2]; ++i) {
          double c1 = xoff + i * xstep;
          double c2 = yoff + j * ystep;
          vec3<double> s1p = cs.to_beam_vector(vec2<double>(c1,c2));
          vec2<double> xyp = panel.get_ray_intersection_px(s1p);
          xyp[0] -= sbox.bbox[0];
          xyp[1] -= sbox.bbox[2];
          xy(j,i) = xyp;
        }
      }

      // Compute the frame numbers of each slice on the grid
      af::shared<double> z(data_.accessor()[0]+1);
      for (std::size_t k = 0; k <= data_.accessor()[0]; ++k) {
        double c3 = zoff + k * zstep;
        double phip = cs.to_rotation_angle_fast(c3);
        z[k] = scan_.get_array_index_from_angle(phip) - sbox.bbox[4];
        if (!(z[k] >= 0 && z[k] <= zs)) {
          return false;
        }
      }

      // Set the mask code
      std::size_t mask_code = Valid | Foreground;

      // Compute the normalized signal in the foreground region
      af::versa< double, af::c_grid<3> > signal(sbox.data.accessor(), 0);
      double total_signal = 0.0;
      for (std::size_t i = 0; i < signal.size(); ++i) {
        std::size_t m = sbox.mask[i];
        if (m & Foreground) {
          if (!(m & Valid)) {
            return false;
          }
          signal[i] = sbox.data[i] - sbox.background[i];
          total_signal += signal[i];
        }
      }
      DIALS_ASSERT(total_signal > 0);
      for (std::size_t i = 0; i < signal.size(); ++i) {
        signal[i] /= total_signal;
      }

      // Get a list of pairs of overlapping polygons
      for (std::size_t j = 0; j < data_.accessor()[1]; ++j) {
        for (std::size_t i = 0; i < data_.accessor()[2]; ++i) {
          vec2<double> xy00 = xy(j,i);
          vec2<double> xy01 = xy(j,i+1);
          vec2<double> xy11 = xy(j+1,i+1);
          vec2<double> xy10 = xy(j+1,i);
          int x0 = std::floor(min4(xy00[0], xy01[0], xy11[0], xy10[0]));
          int x1 = std::ceil (max4(xy00[0], xy01[0], xy11[0], xy10[0]));
          int y0 = std::floor(min4(xy00[1], xy01[1], xy11[1], xy10[1]));
          int y1 = std::ceil (max4(xy00[1], xy01[1], xy11[1], xy10[1]));
          DIALS_ASSERT(x0 < x1);
          DIALS_ASSERT(y0 < y1);
          DIALS_ASSERT(x0 >= 0 && x1 <= xs);
          DIALS_ASSERT(y0 >= 0 && y1 <= ys);
          vert4 p1(xy00, xy01, xy11, xy10);
          reverse_quad_inplace_if_backward(p1);
          for (std::size_t jj = y0; jj < y1; ++jj) {
            for (std::size_t ii = x0; ii < x1; ++ii) {
              vec2<double> xy200(ii,jj);
              vec2<double> xy201(ii,jj+1);
              vec2<double> xy211(ii+1,jj+1);
              vec2<double> xy210(ii+1,jj);
              vert4 p2(xy200, xy201, xy211, xy210);
              reverse_quad_inplace_if_backward(p2);
              vert8 p3 = quad_with_convex_quad(p1, p2);
              double area = simple_area(p3);
              const double EPS = 1e-7;
              DIALS_ASSERT(0.0 <= area && area <= (1.0+EPS));
              if (area > 0) {
                for (std::size_t k = 0; k < data_.accessor()[0]; ++k) {
                  double f00 = std::min(z[k], z[k+1]);
                  double f01 = std::max(z[k], z[k+1]);
                  DIALS_ASSERT(f01 > f00);
                  double fr = f01 - f00;
                  int z0 = std::max((int)0, (int)std::floor(f00));
                  int z1 = std::min((int)zs, (int)std::ceil(f01));
                  DIALS_ASSERT(z0 >= 0 && z1 <= (int)zs);
                  for (int kk = z0; kk < z1; ++kk) {
                    if ((sbox.mask(kk,jj,ii) & mask_code) == mask_code) {
                      std::size_t f10 = kk;
                      std::size_t f11 = kk+1;
                      double fraction = std::min(1.0, (f11 - f10) / fr);
                      DIALS_ASSERT(fraction >= 0.0);
                      double value = fraction * area * signal(kk,jj,ii);
                      data_(k,j,i) += value;
                    }
                  }
                }
              }
            }
          }
        }
      }
      count_++;
      return true;
    }

    /**
     * Finialize the profile.
     */
    void finalize() {
      DIALS_ASSERT(!finalized_);
      finalized_ = true;
    }

    /**
     * @returns The profile.
     */
    profile_type data() const {
      return data_;
    }

    /**
     * @returns The number of reflections used to compute the profile.
     */
    std::size_t count() const {
      return count_;
    }

    profile_type get(
        std::size_t panel_number,
        const vec3<double> &s1,
        double phi,
        int6 bbox) const {

      // Check the input
      DIALS_ASSERT(finalized_);
      DIALS_ASSERT(bbox[1] > bbox[0]);
      DIALS_ASSERT(bbox[3] > bbox[2]);
      DIALS_ASSERT(bbox[5] > bbox[4]);

      // Create the profile array
      std::size_t xs = bbox[1] - bbox[0];
      std::size_t ys = bbox[3] - bbox[2];
      std::size_t zs = bbox[5] - bbox[4];
      profile_type profile(af::c_grid<3>(zs, ys, xs), 1);

      // Compute the grid step and offset
      double xoff = -delta_b_;
      double yoff = -delta_b_;
      double zoff = -delta_m_;
      double xstep = (2.0 * delta_b_) / data_.accessor()[2];
      double ystep = (2.0 * delta_b_) / data_.accessor()[1];
      double zstep = (2.0 * delta_m_) / data_.accessor()[0];

      // Create the coordinate system
      CoordinateSystem cs(m2_, s0_, s1, phi);

      // Get the panel
      const Panel &panel = detector_[panel_number];

      // Compute the detector coordinates of each point on the grid
      af::versa< vec2<double>, af::c_grid<2> > xy(af::c_grid<2>(
            data_.accessor()[1] + 1,
            data_.accessor()[2] + 1));
      for (std::size_t j = 0; j <= data_.accessor()[1]; ++j) {
        for (std::size_t i = 0; i <= data_.accessor()[2]; ++i) {
          double c1 = xoff + i * xstep;
          double c2 = yoff + j * ystep;
          vec3<double> s1p = cs.to_beam_vector(vec2<double>(c1,c2));
          vec2<double> xyp = panel.get_ray_intersection_px(s1p);
          xyp[0] -= bbox[0];
          xyp[1] -= bbox[2];
          xy(j,i) = xyp;
        }
      }

      // Compute the frame numbers of each slice on the grid
      af::shared<double> z(data_.accessor()[0]+1);
      for (std::size_t k = 0; k <= data_.accessor()[0]; ++k) {
        double c3 = zoff + k * zstep;
        double phip = cs.to_rotation_angle_fast(c3);
        z[k] = scan_.get_array_index_from_angle(phip) - bbox[4];
      }

      // Get a list of pairs of overlapping polygons
      for (std::size_t j = 0; j < data_.accessor()[1]; ++j) {
        for (std::size_t i = 0; i < data_.accessor()[2]; ++i) {
          vec2<double> xy00 = xy(j,i);
          vec2<double> xy01 = xy(j,i+1);
          vec2<double> xy11 = xy(j+1,i+1);
          vec2<double> xy10 = xy(j+1,i);
          int x0 = std::floor(min4(xy00[0], xy01[0], xy11[0], xy10[0]));
          int x1 = std::ceil (max4(xy00[0], xy01[0], xy11[0], xy10[0]));
          int y0 = std::floor(min4(xy00[1], xy01[1], xy11[1], xy10[1]));
          int y1 = std::ceil (max4(xy00[1], xy01[1], xy11[1], xy10[1]));
          DIALS_ASSERT(x0 < x1);
          DIALS_ASSERT(y0 < y1);
          DIALS_ASSERT(x0 >= 0 && x1 <= xs);
          DIALS_ASSERT(y0 >= 0 && y1 <= ys);
          vert4 p1(xy00, xy01, xy11, xy10);
          double p1_area = simple_area(p1);
          DIALS_ASSERT(p1_area > 0);
          reverse_quad_inplace_if_backward(p1);
          for (std::size_t jj = y0; jj < y1; ++jj) {
            for (std::size_t ii = x0; ii < x1; ++ii) {
              vec2<double> xy200(ii,jj);
              vec2<double> xy201(ii,jj+1);
              vec2<double> xy211(ii+1,jj+1);
              vec2<double> xy210(ii+1,jj);
              vert4 p2(xy200, xy201, xy211, xy210);
              reverse_quad_inplace_if_backward(p2);
              vert8 p3 = quad_with_convex_quad(p1, p2);
              double area = simple_area(p3);
              area /= p1_area;
              const double EPS = 1e-7;
              DIALS_ASSERT(0.0 <= area && area <= (1.0+EPS));
              if (area > 0) {
                for (std::size_t k = 0; k < data_.accessor()[0]; ++k) {
                  double f00 = std::min(z[k], z[k+1]);
                  double f01 = std::max(z[k], z[k+1]);
                  DIALS_ASSERT(f01 > f00);
                  double fr = f01 - f00;
                  int z0 = std::max((int)0, (int)std::floor(f00));
                  int z1 = std::min((int)zs, (int)std::ceil(f01));
                  DIALS_ASSERT(z0 >= 0 && z1 <= (int)zs);
                  for (int kk = z0; kk < z1; ++kk) {
                    std::size_t f10 = kk;
                    std::size_t f11 = kk+1;
                    double fraction = std::min(1.0, (f11 - f10) / fr);
                    DIALS_ASSERT(fraction >= 0.0);
                    double value = fraction * area * data_(k,j,i);
                    profile(kk,jj,ii) += value;
                  }
                }
              }
            }
          }
        }
      }

      // Return the profile
      return profile;
    }

  private:

    Detector detector_;
    Scan scan_;
    vec3<double> m2_;
    vec3<double> s0_;
    double delta_b_;
    double delta_m_;
    profile_type data_;
    bool finalized_;
    std::size_t count_;
  };


  /**
   * A class to compute a single reference profile per experiment.
   */
  template <typename FloatType>
  class MultiExpSingleProfileLearner {
  public:

    typedef FloatType float_type;
    typedef SingleProfileLearner<FloatType> learner_type;
    typedef typename learner_type::profile_type profile_type;
    typedef typename learner_type::mask_type mask_type;

    /**
     * Compute the reference profiles.
     * @param spec The list of experiment specifications
     * @param data The list of reflections.
     */
    MultiExpSingleProfileLearner(
        af::const_ref<Spec> spec,
        af::reflection_table data,
        std::size_t grid_size) {

      // Check the input
      DIALS_ASSERT(spec.size() > 0);
      DIALS_ASSERT(data.size() > 0);
      DIALS_ASSERT(grid_size > 0);

      // Create the array of learners
      learner_.reserve(spec.size());
      for (std::size_t i = 0; i < spec.size(); ++i) {
        learner_.push_back(learner_type(spec[i],
              int3(grid_size, grid_size, grid_size)));
      }

      // Check the input contains expected fields
      DIALS_ASSERT(data.size() > 0);
      DIALS_ASSERT(data.is_consistent());
      DIALS_ASSERT(data.contains("id"));
      DIALS_ASSERT(data.contains("shoebox"));
      DIALS_ASSERT(data.contains("s1"));
      DIALS_ASSERT(data.contains("xyzcal.mm"));
      DIALS_ASSERT(data.contains("flags"));

      // Get the data we need
      af::const_ref<std::size_t>    id      = data["id"];
      af::const_ref< Shoebox<> >    shoebox = data["shoebox"];
      af::const_ref< vec3<double> > s1      = data["s1"];
      af::const_ref< vec3<double> > xyzcal  = data["xyzcal.mm"];
      af::ref<std::size_t>          flags   = data["flags"];

      // Loop through all the reflections
      for (std::size_t i = 0; i < data.size(); ++i) {
        DIALS_ASSERT(id[i] < learner_.size());
        if (flags[i] & ReferenceSpot) {
          learner_[id[i]].add(shoebox[i], s1[i], xyzcal[i][2]);
        }
      }

      // Finalize the profiles
      for (std::size_t i = 0; i < learner_.size(); ++i) {
        learner_[i].finalize();
      }
    }

    /**
     * Get the reference profile for the reflection.
     * @param id The experiment ID
     * @param panel The panel number
     * @param xyz The calculated xyz pixel coordinate
     * @param bbox The bounding box
     * @returns The reference profile for the reflection.
     */
    profile_type get(
        std::size_t id,
        std::size_t panel,
        const vec3<double> &s1,
        double phi,
        int6 bbox) const {
      DIALS_ASSERT(id < learner_.size());
      return learner_[id].get(panel, s1, phi, bbox);
    }

    /**
     * Get the reference profile.
     * @param id The experiment ID
     * @returns The reference profile.
     */
    profile_type data(std::size_t id) {
      DIALS_ASSERT(id < learner_.size());
      return learner_[id].data();
    }

    /**
     * Get the number of profiles used
     * @param id The experiment ID
     * @returns The number of profiles used.
     */
    std::size_t count(std::size_t id) {
      DIALS_ASSERT(id < learner_.size());
      return learner_[id].count();
    }

    /**
     * Return the number of references.
     */
    std::size_t size() const {
      return learner_.size();
    }

  private:

    std::vector<learner_type> learner_;
  };


  /**
   * A class to perform image space profile fitting.
   */
  class ImageSpaceProfileFitting {
  public:

    typedef Shoebox<>::float_type float_type;
    typedef MultiExpSingleProfileLearner<float_type> reference_learner_type;
    typedef reference_learner_type::profile_type profile_type;
    typedef af::versa< bool, af::c_grid<3> > mask_type;
    typedef ProfileFitting<float_type> fitting_type;

    ImageSpaceProfileFitting(std::size_t grid_size)
        : grid_size_(grid_size) {
      DIALS_ASSERT(grid_size > 0);
    }

    /**
     * Add some experiment data
     * @param spec The experiment specification
     */
    void add(const Spec &spec) {
      spec_.push_back(spec);
    }

    /**
     * Perform the profile fitting on the input data
     * @param data The reflection data
     */
    reference_learner_type
    execute(af::reflection_table data) const {

      // Check we have experimental setup data
      DIALS_ASSERT(spec_.size() > 0);

      // Check the input contains expected fields
      DIALS_ASSERT(data.size() > 0);
      DIALS_ASSERT(data.is_consistent());
      DIALS_ASSERT(data.contains("id"));
      DIALS_ASSERT(data.contains("shoebox"));
      DIALS_ASSERT(data.contains("s1"));
      DIALS_ASSERT(data.contains("xyzcal.mm"));
      DIALS_ASSERT(data.contains("flags"));

      // Compute the reference profile
      reference_learner_type reference(spec_.const_ref(), data, grid_size_);

      // Get the data we need
      af::const_ref<std::size_t>    id      = data["id"];
      af::const_ref< Shoebox<> >    shoebox = data["shoebox"];
      af::const_ref< vec3<double> > s1      = data["s1"];
      af::const_ref< vec3<double> > xyzcal  = data["xyzcal.mm"];
      af::ref<std::size_t>          flags   = data["flags"];

      // Get the new columns to set
      af::ref<double> intensity   = data["intensity.prf.value"];
      af::ref<double> variance    = data["intensity.prf.variance"];
      af::ref<double> correlation = data["profile_correlation"];

      // Do the profile fitting for all reflections
      for (std::size_t i = 0; i < data.size(); ++i) {

        // Initialise some stuff
        intensity[i] = 0;
        variance[i] = 0;
        correlation[i] = 0;
        flags[i] &= ~IntegratedPrf;

        // Integrate
        if (!(flags[i] & DontIntegrate)) {

          // Get the shoebox
          const Shoebox<> &sbox = shoebox[i];
          DIALS_ASSERT(sbox.is_consistent());

          // Get the profile for a given reflection
          profile_type profile = reference.get(
              id[i],
              sbox.panel,
              s1[i],
              xyzcal[i][2],
              sbox.bbox);

          // Compute the integration mask
          mask_type mask(sbox.mask.accessor(), false);
          std::size_t mask_code = Valid | Foreground;
          std::size_t mask_count = 0;
          DIALS_ASSERT(profile.accessor().all_eq(sbox.mask.accessor()));
          for (std::size_t j = 0; j < sbox.mask.size(); ++j) {
            if ((sbox.mask[j] & mask_code) == mask_code && profile[j] >= 0.0) {
              mask[j] = true;
              mask_count++;
            }
          }

          // Perform the profile fit
          if (mask_count > 0) {
            fitting_type fit(
                profile.const_ref(),
                mask.const_ref(),
                sbox.data.const_ref(),
                sbox.background.const_ref());

            // Set the data in the reflection
            intensity[i]   = fit.intensity();
            variance[i]    = fit.variance();
            correlation[i] = fit.correlation();

            // Set the integrated flag
            flags[i] |= IntegratedPrf;
          }
        }
      }

      // Return the reference learner
      return reference;
    }

  private:

    af::shared<Spec> spec_;
    std::size_t grid_size_;
  };

}} // namespace dials::algorithms


#endif // DIALS_ALGORITHMS_INTEGRATION_FIT_IMAGE_FIT_H
