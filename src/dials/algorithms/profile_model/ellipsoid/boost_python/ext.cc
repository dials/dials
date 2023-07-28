/*
 * ext.cc
 *
 *  Copyright (C) 2018 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <scitbx/constants.h>
#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include <scitbx/mat2.h>
#include <scitbx/mat3.h>
#include <scitbx/sym_mat2.h>
#include <scitbx/math/r3_rotation.h>
#include <scitbx/matrix/multiply.h>
#include <scitbx/matrix/eigensystem.h>
#include <dxtbx/model/experiment.h>
#include <dials/model/data/shoebox.h>
#include <dials/algorithms/profile_model/gaussian_rs/coordinate_system.h>
#include <dials/algorithms/shoebox/mask_code.h>
#include <dials/array_family/reflection_table.h>
#include <dials/error.h>

using namespace boost::python;

namespace dials { namespace algorithms { namespace boost_python {

  using dials::algorithms::profile_model::gaussian_rs::CoordinateSystem2d;
  using dials::model::Background;
  using dials::model::Foreground;
  using dials::model::Overlapped;
  using dials::model::Shoebox;
  using dials::model::Valid;
  using dxtbx::model::BeamBase;
  using dxtbx::model::Detector;
  using dxtbx::model::Experiment;
  using dxtbx::model::Panel;
  using scitbx::mat2;
  using scitbx::mat3;
  using scitbx::sym_mat2;
  using scitbx::vec2;
  using scitbx::vec3;
  using scitbx::af::int6;
  using scitbx::matrix::multiply_transpose;
  using scitbx::matrix::transpose_multiply;

  namespace detail {

    /**
     * Helper function to do matrix multiplication
     */
    double AT_B_A(vec2<double> A, mat2<double> B) {
      vec2<double> ATB = A * B;
      return ATB * A;
    }

    /**
     * Helper function to do matrix multiplication
     */
    double AT_B_A(vec3<double> A, mat3<double> B) {
      vec3<double> ATB = A * B;
      return ATB * A;
    }

  }  // namespace detail

  /**
   * Function to given chisq quantile value
   * @param k The degrees of freedom
   * @param p The probability
   */
  double chisq_pdf(int k, double x) {
    DIALS_ASSERT(k > 0);
    DIALS_ASSERT(x >= 0);
    boost::math::chi_squared_distribution<> dist(k);
    return boost::math::pdf(dist, x);
  }

  /**
   * Function to given chisq quantile value
   * @param k The degrees of freedom
   * @param p The probability
   */
  double chisq_quantile(int k, double p) {
    DIALS_ASSERT(k > 0);
    DIALS_ASSERT(p >= 0 && p <= 1);
    boost::math::chi_squared_distribution<> dist(k);
    return boost::math::quantile(dist, p);
  }

  /**
   * Perform the change of basis from lab coordinates to the kabsch coordinate
   * system
   * @param s0 The incident beam vector
   * @param s2 The diffracted beam vector
   * @return The change of basis matrix
   */
  mat3<double> compute_change_of_basis_operation(vec3<double> s0, vec3<double> s2) {
    const double TINY = 1e-7;
    DIALS_ASSERT((s2 - s0).length() > TINY);
    vec3<double> e1 = s2.cross(s0).normalize();
    vec3<double> e2 = s2.cross(e1).normalize();
    vec3<double> e3 = s2.normalize();
    mat3<double> R(e1[0], e1[1], e1[2], e2[0], e2[1], e2[2], e3[0], e3[1], e3[2]);
    return R;
  }

  boost::python::tuple reflection_statistics(const Panel panel,
                                             const vec3<double> xyzobs,
                                             const double s0_length,
                                             const vec3<double> s0,
                                             const Shoebox<> sbox) {
    typedef Shoebox<>::float_type float_type;
    vec2<double> p1(xyzobs[0], xyzobs[1]);
    vec3<double> sp = panel.get_pixel_lab_coord(p1).normalize() * s0_length;
    mat3<double> R = compute_change_of_basis_operation(s0, sp);

    const af::versa<int, af::c_grid<3>> mask = sbox.mask;
    const af::const_ref<float_type, af::c_grid<3>> data = sbox.data.const_ref();
    const af::versa<float_type, af::c_grid<3>> bgrd = sbox.background;
    int i0 = sbox.bbox[0];
    int j0 = sbox.bbox[2];
    int n1 = data.accessor()[1];
    int n2 = data.accessor()[2];

    af::versa<double, af::c_grid<3>> X(af::c_grid<3>(n1, n2, 2));
    af::versa<double, af::c_grid<2>> C(af::c_grid<2>(n1, n2));
    double ctot = 0;
    for (int j = 0; j < n1; ++j) {
      for (int i = 0; i < n2; ++i) {
        double c = data(0, j, i) - bgrd(0, j, i);
        if (c > 0) {
          if ((mask(0, j, i) & (1 | 4)) == (1 | 4)) {
            ctot += c;
            int ii = i + i0;
            int jj = j + j0;
            vec2<double> sc(ii + 0.5, jj + 0.5);
            vec3<double> s = panel.get_pixel_lab_coord(sc).normalize() * s0_length;
            vec3<double> e = R * s;
            X(j, i, 0) = e[0];
            X(j, i, 1) = e[1];
            C(j, i) = c;
          }
        }
      }
    }
    // check ctot > 0

    // now sum to get Xbar
    vec2<double> xbar(0, 0);
    for (int j = 0; j < n1; ++j) {
      for (int i = 0; i < n2; ++i) {
        xbar[0] += X(j, i, 0) * C(j, i);
        xbar[1] += X(j, i, 1) * C(j, i);
      }
    }
    xbar[0] = xbar[0] / ctot;
    xbar[1] = xbar[1] / ctot;

    mat2<double> Sobs(0, 0, 0, 0);
    for (int j = 0; j < n1; ++j) {
      for (int i = 0; i < n2; ++i) {
        double c_i = C(j, i);
        vec2<double> x_i(X(j, i, 0) - xbar[0], X(j, i, 1) - xbar[1]);
        Sobs[0] += pow(x_i[0], 2) * c_i;
        Sobs[1] += x_i[0] * x_i[1] * c_i;
        Sobs[2] += x_i[0] * x_i[1] * c_i;
        Sobs[3] += pow(x_i[1], 2) * c_i;
      }
    }
    Sobs[0] = Sobs[0] / ctot;
    Sobs[1] = Sobs[1] / ctot;
    Sobs[2] = Sobs[2] / ctot;
    Sobs[3] = Sobs[3] / ctot;

    return boost::python::make_tuple(sp, ctot, xbar, Sobs);
  }

  boost::python::list calc_s1_s2(const af::shared<cctbx::miller::index<>> miller_idx,
                                 const mat3<double> A,
                                 const vec3<double> s0) {
    double len_s0 = s0.length();
    af::shared<vec3<double>> s1(miller_idx.size());
    af::shared<vec3<double>> s2(miller_idx.size());
    for (std::size_t i = 0; i < miller_idx.size(); ++i) {
      vec3<double> r = A * miller_idx[i];
      vec3<double> s2_i = r + s0;
      s2[i] = s2_i;
      s1[i] = len_s0 * s2_i.normalize();
    }
    boost::python::list result;
    result.append(s1);
    result.append(s2);
    return result;
  }

  /**
   * Perform the change of basis from reciprocal space to a coordinate system of
   * orientated with the reciprocal lattice vector to the reflection
   * @param s0 The incident beam vector
   * @param r The reciprocal lattice vector
   * @return The change of basis matrix
   */
  mat3<double> compute_change_of_basis_operation2(vec3<double> s0, vec3<double> r) {
    vec3<double> e1 = r.cross(s0).normalize();
    vec3<double> e2 = r.cross(e1).normalize();
    vec3<double> e3 = r.normalize();
    mat3<double> R(e1[0], e1[1], e1[2], e2[0], e2[1], e2[2], e3[0], e3[1], e3[2]);
    return R;
  }

  /**
   * A class to predict the reflections
   */
  class PredictorBase {
  public:
    /**
     * Initialise the predictor
     * @param experiment The experiment
     * @param probability The probability
     */
    PredictorBase(Experiment experiment, double probability)
        : experiment_(experiment), probability_(probability) {
      DIALS_ASSERT(probability > 0 && probability < 1);
    }

    /**
     * Predict the reflections
     * @param h The list of miller indices
     * @returns The reflection table
     */
    af::reflection_table predict(af::const_ref<cctbx::miller::index<>> h) const {
      // Get the beam and detector
      const double TINY = 1e-7;
      DIALS_ASSERT(experiment_.get_beam() != NULL);
      DIALS_ASSERT(experiment_.get_crystal() != NULL);
      DIALS_ASSERT(experiment_.get_detector() != NULL);
      Detector detector = *experiment_.get_detector();

      // Compute quantile
      double quantile = chisq_quantile(3, probability_);

      // Get stuff from experiment
      mat3<double> A = experiment_.get_crystal()->get_A();
      vec3<double> s0 = experiment_.get_beam()->get_s0();

      // Initialise some arrays
      af::shared<cctbx::miller::index<>> miller_indices;
      af::shared<bool> entering;
      af::shared<vec3<double>> s1_list;
      af::shared<vec3<double>> s2_list;
      af::shared<vec3<double>> xyzcalpx;
      af::shared<vec3<double>> xyzcalmm;
      af::shared<std::size_t> panel_list;
      af::shared<int> experiment_id;

      // Loop through the input miller indices
      for (std::size_t i = 0; i < h.size(); ++i) {
        // Compute the point and the distance from the Ewald sphere
        vec3<double> r = A * h[i];
        vec3<double> s2 = s0 + r;
        vec3<double> s3 = s2.normalize() * s0.length();

        // Invert the matrix
        mat3<double> sigma = get_sigma(s0, r);
        mat3<double> sigma_inv = sigma.inverse();

        // Compute distance
        double d = detail::AT_B_A(s3 - s2, sigma_inv);

        // If it is close enough then predict stuff
        if (d < quantile) {
          // Compute the rotation of the reflection
          mat3<double> R = compute_change_of_basis_operation(s0, s2);

          // Rotate the covariance matrix and s2 vector
          mat3<double> S = R * sigma * R.transpose();
          vec3<double> mu = R * s2;
          vec3<double> zaxis(0, 0, 1);
          DIALS_ASSERT(std::abs((mu.normalize() * zaxis) - 1) < TINY);

          // Partition the covariance matrix
          mat2<double> S11(S[0], S[1], S[3], S[4]);
          vec2<double> S12(S[2], S[5]);
          vec2<double> S21(S[6], S[7]);
          double S22 = S[8];

          // Partition the mean vector
          vec2<double> mu1(mu[0], mu[1]);
          double mu2 = mu[2];

          // Compute epsilon the distance to the Ewald sphere
          double epsilon = s0.length() - mu2;

          // Compute the mean of the conditional distribution
          DIALS_ASSERT(S22 > 0);
          double S22_inv = 1.0 / S22;
          vec2<double> mubar = mu1 + S12 * S22_inv * epsilon;

          // Compute the diffracted beam vector
          vec3<double> v(mubar[0], mubar[1], s0.length());
          vec3<double> s1 = R.transpose() * (v.normalize() * s0.length());
          try {
            // Do the panel ray intersection
            Detector::coord_type impact = detector.get_ray_intersection(s1);
            std::size_t panel_id = impact.first;
            Panel panel = detector[panel_id];
            vec2<double> xymm = panel.get_ray_intersection(s1);
            vec2<double> xypx = panel.millimeter_to_pixel(xymm);

            // Append the stuff to arrays
            experiment_id.push_back(0);
            miller_indices.push_back(h[i]);
            entering.push_back(false);
            panel_list.push_back(panel_id);
            s1_list.push_back(s1);
            s2_list.push_back(s2);
            xyzcalpx.push_back(vec3<double>(xypx[0], xypx[1], 0));
            xyzcalmm.push_back(vec3<double>(xymm[0], xymm[1], 0));
          } catch (dxtbx::error) {
            continue;
          }
        }
      }

      // Construct the reflection table
      af::reflection_table reflections(miller_indices.size());
      reflections["miller_index"] = miller_indices;
      reflections["entering"] = entering;
      reflections["s1"] = s1_list;
      reflections["s2"] = s2_list;
      reflections["xyzcal.px"] = xyzcalpx;
      reflections["xyzcal.mm"] = xyzcalmm;
      reflections["panel"] = panel_list;
      reflections["id"] = experiment_id;
      return reflections;
    }

  protected:
    /**
     * Get the sigma for a reflection
     */
    virtual mat3<double> get_sigma(vec3<double> s0, vec3<double> r) const {
      throw DIALS_ERROR("Overload!");
    }

    Experiment experiment_;
    double probability_;
  };

  /**
   * The predictor for simple profile models
   */
  class PredictorSimple : public PredictorBase {
  public:
    /**
     * Initialise the class
     * @param experiment The experiment
     * @param sigma The covariance matrix
     */
    PredictorSimple(Experiment experiment, mat3<double> sigma, double probability)
        : PredictorBase(experiment, probability), sigma_(sigma) {}

  protected:
    /**
     * Get the sigma for a reflection
     */
    virtual mat3<double> get_sigma(vec3<double> s0, vec3<double> r) const {
      return sigma_;
    }

    mat3<double> sigma_;
  };

  /**
   * The predictor for angular profile models
   */
  class PredictorAngular : public PredictorBase {
  public:
    /**
     * Initialise the class
     * @param experiment The experiment
     * @param sigma The covariance matrix
     */
    PredictorAngular(Experiment experiment,
                     mat3<double> sigma,
                     double probability,
                     mat3<double> sigma_A)
        : PredictorBase(experiment, probability), sigma_(sigma), sigma_A_(sigma_A) {}

  protected:
    /**
     * Get the sigma for a reflection
     */
    virtual mat3<double> get_sigma(vec3<double> s0, vec3<double> r) const {
      const double TINY = 1e-7;
      mat3<double> Q = compute_change_of_basis_operation2(s0, r);
      vec3<double> zaxis(0, 0, 1);
      DIALS_ASSERT(std::abs(((Q * r.normalize()) * zaxis) - 1) < TINY);
      double r_scale = std::pow(r.length(), 2);
      mat3<double> A(r_scale, 0, 0, 0, r_scale, 0, 0, 0, 0.0);
      return (Q.transpose() * A * sigma_A_ * Q) + sigma_;
    }

    mat3<double> sigma_;
    mat3<double> sigma_A_;
  };

  /**
   * A class to compute the bounding box
   */
  class BBoxCalculatorBase {
  public:
    /**
     * Initialise the class
     * @param experiment The experiment
     * @param sigma The covariance matrix
     */
    BBoxCalculatorBase(Experiment experiment, double probability, int border)
        : experiment_(experiment), probability_(probability), border_(border) {
      DIALS_ASSERT(border > 0);
      DIALS_ASSERT(probability < 1.0);
      DIALS_ASSERT(probability > 0);
    }

    /**
     * Compute the bounding box
     * @param reflections The reflection table
     */
    void compute(af::reflection_table reflections) const {
      // Get some array from the reflection table
      af::const_ref<vec3<double>> s1 = reflections["s1"];
      af::const_ref<vec3<double>> s2 = reflections["s2"];
      af::const_ref<std::size_t> panel = reflections["panel"];
      af::ref<int6> bbox = reflections["bbox"];

      // Compute quantile
      double D = chisq_quantile(2, probability_);

      // Compute mask for all reflections
      for (std::size_t i = 0; i < reflections.size(); ++i) {
        bbox[i] = compute_single(s1[i], s2[i], panel[i], D);
      }
    }

  protected:
    /**
     * Get the sigma for a reflection
     */
    virtual mat3<double> get_sigma(vec3<double> s0, vec3<double> r) const {
      throw DIALS_ERROR("Overload!");
    }

    /**
     * Compute the bounding box for a single reflection
     * @param s1 The diffracted beam vector touching the Ewald sphere
     * @param s2 The diffracted beam vector to the centre of the rlp
     * @param panel_id The panel id
     * @param D the distance
     */
    int6 compute_single(vec3<double> s1,
                        vec3<double> s2,
                        std::size_t panel_id,
                        double D) const {
      const double TINY = 1e-7;
      DIALS_ASSERT(s1.length() > 0);
      DIALS_ASSERT(s2.length() > 0);
      DIALS_ASSERT(D > 0);

      // Get the beam and detector
      DIALS_ASSERT(experiment_.get_beam() != NULL);
      DIALS_ASSERT(experiment_.get_detector() != NULL);
      Detector detector = *experiment_.get_detector();

      // The the indicident beam vector
      vec3<double> s0 = experiment_.get_beam()->get_s0();
      double s0_length = s0.length();
      DIALS_ASSERT(std::abs(s0_length - s1.length()) < TINY);

      vec3<double> r = s2 - s0;

      // Compute the change of basis for the reflection
      mat3<double> R = compute_change_of_basis_operation(s0, s2);

      // Rotate the covariance matrix and s2 vector
      mat3<double> S = R * get_sigma(s0, r) * R.transpose();
      vec3<double> mu = R * s2;
      vec3<double> zaxis(0, 0, 1);
      DIALS_ASSERT(std::abs((mu.normalize() * zaxis) - 1) < TINY);

      // Partition the covariance matrix
      mat2<double> S11(S[0], S[1], S[3], S[4]);
      vec2<double> S12(S[2], S[5]);
      vec2<double> S21(S[6], S[7]);
      double S22 = S[8];

      // Partition the mean vector
      vec2<double> mu1(mu[0], mu[1]);
      double mu2 = mu[2];

      // Compute epsilon the distance to the Ewald sphere
      double epsilon = s0.length() - mu2;

      // Compute the mean of the conditional distribution
      DIALS_ASSERT(S22 > 0);
      double S22_inv = 1.0 / S22;
      vec2<double> mubar = mu1 + S12 * S22_inv * epsilon;

      // Compute the covariance of the conditional distribution
      mat2<double> S12_S21;
      multiply_transpose(&S12[0], &S21[0], 2, 1, 2, &S12_S21[0]);
      mat2<double> Sbar = S11 - S12_S21 * S22_inv;

      // Get the panel model
      Panel panel = detector[panel_id];

      // Compute the the min/max bounding box on the ellipse
      double A1 = Sbar[0];
      double A2 = Sbar[3];
      DIALS_ASSERT(A1 >= 0 && A2 >= 0);
      double delta1 = std::sqrt(D * A1);
      double delta2 = std::sqrt(D * A2);

      // The corner points in conditional space
      vec2<double> p1 = mubar + vec2<double>(-delta1, -delta2);
      vec2<double> p2 = mubar + vec2<double>(-delta1, +delta2);
      vec2<double> p3 = mubar + vec2<double>(+delta1, -delta2);
      vec2<double> p4 = mubar + vec2<double>(+delta1, +delta2);

      // The corner points in lab space
      vec3<double> sp1 = R.transpose() * vec3<double>(p1[0], p1[1], s0.length());
      vec3<double> sp2 = R.transpose() * vec3<double>(p2[0], p2[1], s0.length());
      vec3<double> sp3 = R.transpose() * vec3<double>(p3[0], p3[1], s0.length());
      vec3<double> sp4 = R.transpose() * vec3<double>(p4[0], p4[1], s0.length());

      // The xy coordinates on the detector
      vec2<double> xy1 = panel.get_ray_intersection_px(sp1);
      vec2<double> xy2 = panel.get_ray_intersection_px(sp2);
      vec2<double> xy3 = panel.get_ray_intersection_px(sp3);
      vec2<double> xy4 = panel.get_ray_intersection_px(sp4);

      // Get the min and max x and y coords
      double xmin = std::min(std::min(xy1[0], xy2[0]), std::min(xy3[0], xy4[0]));
      double ymin = std::min(std::min(xy1[1], xy2[1]), std::min(xy3[1], xy4[1]));
      double xmax = std::max(std::max(xy1[0], xy2[0]), std::max(xy3[0], xy4[0]));
      double ymax = std::max(std::max(xy1[1], xy2[1]), std::max(xy3[1], xy4[1]));

      // Create bounding box
      int x0 = ((int)std::floor(xmin)) - border_;
      int y0 = ((int)std::floor(ymin)) - border_;
      int x1 = ((int)std::ceil(xmax)) + border_;
      int y1 = ((int)std::ceil(ymax)) + border_;
      DIALS_ASSERT(x1 > x0);
      DIALS_ASSERT(y1 > y0);
      return int6(x0, x1, y0, y1, 0, 1);
    }

    Experiment experiment_;
    double probability_;
    int border_;
  };

  /**
   * The bbox calculator for simple profile models
   */
  class BBoxCalculatorSimple : public BBoxCalculatorBase {
  public:
    /**
     * Initialise the class
     * @param experiment The experiment
     * @param sigma The covariance matrix
     */
    BBoxCalculatorSimple(Experiment experiment,
                         mat3<double> sigma,
                         double probability,
                         int border)
        : BBoxCalculatorBase(experiment, probability, border), sigma_(sigma) {}

  protected:
    /**
     * Get the sigma for a reflection
     */
    virtual mat3<double> get_sigma(vec3<double> s0, vec3<double> r) const {
      return sigma_;
    }

    mat3<double> sigma_;
  };

  /**
   * The bbox calculator for angular profile models
   */
  class BBoxCalculatorAngular : public BBoxCalculatorBase {
  public:
    /**
     * Initialise the class
     * @param experiment The experiment
     * @param sigma The covariance matrix
     */
    BBoxCalculatorAngular(Experiment experiment,
                          mat3<double> sigma,
                          double probability,
                          int border,
                          mat3<double> sigma_A)
        : BBoxCalculatorBase(experiment, probability, border),
          sigma_(sigma),
          sigma_A_(sigma_A) {}

  protected:
    /**
     * Get the sigma for a reflection
     */
    virtual mat3<double> get_sigma(vec3<double> s0, vec3<double> r) const {
      const double TINY = 1e-7;
      mat3<double> Q = compute_change_of_basis_operation2(s0, r);
      vec3<double> zaxis(0, 0, 1);
      DIALS_ASSERT(std::abs(((Q * r.normalize()) * zaxis) - 1) < TINY);
      double r_scale = std::pow(r.length(), 2);
      mat3<double> A(r_scale, 0, 0, 0, r_scale, 0, 0, 0, 0.0);
      return (Q.transpose() * A * sigma_A_ * Q) + sigma_;
    }

    mat3<double> sigma_;
    mat3<double> sigma_A_;
  };

  /**
   * A class to compute the reflection mask
   */
  class MaskCalculatorBase {
  public:
    /**
     * Initialise the class
     * @param experiment The experiment
     * @param sigma The covariance matrix
     */
    MaskCalculatorBase(Experiment experiment, double probability)
        : experiment_(experiment), probability_(probability) {
      DIALS_ASSERT(probability < 1.0);
      DIALS_ASSERT(probability > 0);
    }

    /**
     * Compute the reflection mask
     * @param reflections The reflection table
     */
    void compute(af::reflection_table reflections) const {
      // Get some array from the reflection table
      af::const_ref<vec3<double>> s1 = reflections["s1"];
      af::const_ref<vec3<double>> s2 = reflections["s2"];
      af::ref<Shoebox<>> sbox = reflections["shoebox"];

      // Compute quantile
      double D = chisq_quantile(2, probability_);

      // Compute mask for all reflections
      for (std::size_t i = 0; i < reflections.size(); ++i) {
        compute_single(s1[i], s2[i], sbox[i], D);
      }
    }

  protected:
    /**
     * Get the sigma for a reflection
     */
    virtual mat3<double> get_sigma(vec3<double> s0, vec3<double> r) const {
      throw DIALS_ERROR("Overload!");
    }

    /**
     * Compute the bounding box for a single reflection
     * @param s1 The diffracted beam vector touching the Ewald sphere
     * @param s2 The diffracted beam vector to the centre of the rlp
     * @param sbox The shoebox
     * @param D the distance
     */
    void compute_single(vec3<double> s1,
                        vec3<double> s2,
                        Shoebox<> sbox,
                        double D) const {
      const double TINY = 1e-7;
      DIALS_ASSERT(s1.length() > 0);
      DIALS_ASSERT(s2.length() > 0);
      DIALS_ASSERT(sbox.is_consistent());
      DIALS_ASSERT(D > 0);

      // Get the beam and detector
      DIALS_ASSERT(experiment_.get_beam() != NULL);
      DIALS_ASSERT(experiment_.get_detector() != NULL);
      Detector detector = *experiment_.get_detector();

      // The the indicident beam vector
      vec3<double> s0 = experiment_.get_beam()->get_s0();
      double s0_length = s0.length();
      DIALS_ASSERT(std::abs(s0_length - s1.length()) < TINY);

      vec3<double> r = s2 - s0;

      // Compute the change of basis for the reflection
      mat3<double> R = compute_change_of_basis_operation(s0, s2);

      // Rotate the covariance matrix and s2 vector
      mat3<double> S = R * get_sigma(s0, r) * R.transpose();
      vec3<double> mu = R * s2;
      vec3<double> zaxis(0, 0, 1);
      DIALS_ASSERT(std::abs((mu.normalize() * zaxis) - 1) < TINY);

      // Partition the covariance matrix
      mat2<double> S11(S[0], S[1], S[3], S[4]);
      vec2<double> S12(S[2], S[5]);
      vec2<double> S21(S[6], S[7]);
      double S22 = S[8];

      // Partition the mean vector
      vec2<double> mu1(mu[0], mu[1]);
      double mu2 = mu[2];

      // Compute epsilon the distance to the Ewald sphere
      double epsilon = s0.length() - mu2;

      // Compute the mean of the conditional distribution
      DIALS_ASSERT(S22 > 0);
      double S22_inv = 1.0 / S22;
      vec2<double> mubar = mu1 + S12 * S22_inv * epsilon;

      // Compute the covariance of the conditional distribution
      mat2<double> S12_S21;
      multiply_transpose(&S12[0], &S21[0], 2, 1, 2, &S12_S21[0]);
      mat2<double> Sbar = S11 - S12_S21 * S22_inv;
      mat2<double> Sbar_inv = Sbar.inverse();

      // Get the mask array
      af::ref<int, af::c_grid<3>> mask = sbox.mask.ref();

      // Get the bounding box
      int x0 = sbox.bbox[0];
      int x1 = sbox.bbox[1];
      int y0 = sbox.bbox[2];
      int y1 = sbox.bbox[3];
      int z0 = sbox.bbox[4];
      int z1 = sbox.bbox[5];
      DIALS_ASSERT(x0 < x1);
      DIALS_ASSERT(y0 < y1);
      DIALS_ASSERT(z0 < z1);
      DIALS_ASSERT(z1 - z0 == 1);

      // Create the coordinate system
      CoordinateSystem2d cs(s0, s2);

      // Get the panel model
      Panel panel = detector[sbox.panel];

      // Set the mask value for each pixel
      DIALS_ASSERT(mask.accessor()[0] == 1);
      for (std::size_t j = 0; j < mask.accessor()[1]; ++j) {
        for (std::size_t i = 0; i < mask.accessor()[2]; ++i) {
          int ii = x0 + ((int)i);
          int jj = y0 + ((int)j);

          // The pixel coordinates of the corners
          vec2<double> p1(ii, jj);
          vec2<double> p2(ii + 1, jj);
          vec2<double> p3(ii, jj + 1);
          vec2<double> p4(ii + 1, jj + 1);

          // The lab coordinates of the pixel corners
          vec3<double> sp1 = panel.get_pixel_lab_coord(p1).normalize() * s0_length;
          vec3<double> sp2 = panel.get_pixel_lab_coord(p2).normalize() * s0_length;
          vec3<double> sp3 = panel.get_pixel_lab_coord(p3).normalize() * s0_length;
          vec3<double> sp4 = panel.get_pixel_lab_coord(p4).normalize() * s0_length;

          // The coordinates in kabsch space
          vec2<double> x1 = cs.from_beam_vector(sp1);
          vec2<double> x2 = cs.from_beam_vector(sp2);
          vec2<double> x3 = cs.from_beam_vector(sp3);
          vec2<double> x4 = cs.from_beam_vector(sp4);

          // The distance from the mean
          double d1 = detail::AT_B_A(x1 - mubar, Sbar_inv);
          double d2 = detail::AT_B_A(x2 - mubar, Sbar_inv);
          double d3 = detail::AT_B_A(x3 - mubar, Sbar_inv);
          double d4 = detail::AT_B_A(x4 - mubar, Sbar_inv);

          // The minimum distance
          if (std::min(std::min(d1, d2), std::min(d3, d4)) < D) {
            mask(0, j, i) |= Foreground;
          } else {
            mask(0, j, i) |= Background;
          }
        }
      }
    }

    Experiment experiment_;
    double probability_;
  };

  /**
   * The mask calculator for simple profile models
   */
  class MaskCalculatorSimple : public MaskCalculatorBase {
  public:
    /**
     * Initialise the class
     * @param experiment The experiment
     * @param sigma The covariance matrix
     */
    MaskCalculatorSimple(Experiment experiment, mat3<double> sigma, double probability)
        : MaskCalculatorBase(experiment, probability), sigma_(sigma) {}

  protected:
    /**
     * Get the sigma for a reflection
     */
    virtual mat3<double> get_sigma(vec3<double> s0, vec3<double> r) const {
      return sigma_;
    }

    mat3<double> sigma_;
  };

  /**
   * The mask calculator for angular profile models
   */
  class MaskCalculatorAngular : public MaskCalculatorBase {
  public:
    /**
     * Initialise the class
     * @param experiment The experiment
     * @param sigma The covariance matrix
     */
    MaskCalculatorAngular(Experiment experiment,
                          mat3<double> sigma,
                          double probability,
                          mat3<double> sigma_A)
        : MaskCalculatorBase(experiment, probability),
          sigma_(sigma),
          sigma_A_(sigma_A) {}

  protected:
    /**
     * Get the sigma for a reflection
     */
    virtual mat3<double> get_sigma(vec3<double> s0, vec3<double> r) const {
      const double TINY = 1e-7;
      mat3<double> Q = compute_change_of_basis_operation2(s0, r);
      vec3<double> zaxis(0, 0, 1);
      DIALS_ASSERT(std::abs(((Q * r.normalize()) * zaxis) - 1) < TINY);
      double r_scale = std::pow(r.length(), 2);
      mat3<double> A(r_scale, 0, 0, 0, r_scale, 0, 0, 0, 0.0);
      return (Q.transpose() * A * sigma_A_ * Q) + sigma_;
    }

    mat3<double> sigma_;
    mat3<double> sigma_A_;
  };

  vec2<double> rse(const mat3<double> &R,
                   const vec2<double> &mbar,
                   const vec2<double> &xobs,
                   const float norm_s0,
                   const Detector &detector) {
    // double norm_s0 = 1.0;
    vec3<double> s1;
    vec3<double> s3;
    s1[0] = (R[0] * mbar[0]) + (R[3] * mbar[1]) + (R[6] * norm_s0);
    s1[1] = (R[1] * mbar[0]) + (R[4] * mbar[1]) + (R[7] * norm_s0);
    s1[2] = (R[2] * mbar[0]) + (R[5] * mbar[1]) + (R[8] * norm_s0);
    s3[0] = (R[0] * xobs[0]) + (R[3] * xobs[1]) + (R[6] * norm_s0);
    s3[1] = (R[1] * xobs[0]) + (R[4] * xobs[1]) + (R[7] * norm_s0);
    s3[2] = (R[2] * xobs[0]) + (R[5] * xobs[1]) + (R[8] * norm_s0);
    vec2<double> xyzcal = detector[0].get_ray_intersection_px(s1);
    vec2<double> xyzobs = detector[0].get_ray_intersection_px(s3);
    double rx2 = pow(xyzcal[0] - xyzobs[0], 2);
    double ry2 = pow(xyzcal[1] - xyzobs[1], 2);
    vec2<double> res{rx2, ry2};
    return res;
  }

  BOOST_PYTHON_MODULE(dials_algorithms_profile_model_ellipsoid_ext) {
    def("chisq_quantile", &chisq_quantile);
    def("chisq_pdf", &chisq_pdf);
    def("reflection_statistics", &reflection_statistics);
    def("calc_s1_s2", &calc_s1_s2);
    def("rse", &rse);

    class_<PredictorBase>("PredictorBase", no_init)
      .def("predict", &PredictorBase::predict);

    class_<PredictorSimple, bases<PredictorBase>>("PredictorSimple", no_init)
      .def(init<Experiment, mat3<double>, double>());

    class_<PredictorAngular, bases<PredictorBase>>("PredictorAngular", no_init)
      .def(init<Experiment, mat3<double>, double, mat3<double>>());

    class_<BBoxCalculatorBase>("BBoxCalculatorBase", no_init)
      .def("compute", &BBoxCalculatorBase::compute);

    class_<BBoxCalculatorSimple, bases<BBoxCalculatorBase>>("BBoxCalculatorSimple",
                                                            no_init)
      .def(init<Experiment, mat3<double>, double, int>());

    class_<BBoxCalculatorAngular, bases<BBoxCalculatorBase>>("BBoxCalculatorAngular",
                                                             no_init)
      .def(init<Experiment, mat3<double>, double, int, mat3<double>>());

    class_<MaskCalculatorBase>("MaskCalculatorBase", no_init)
      .def("compute", &MaskCalculatorBase::compute);

    class_<MaskCalculatorSimple, bases<MaskCalculatorBase>>("MaskCalculatorSimple",
                                                            no_init)
      .def(init<Experiment, mat3<double>, double>());

    class_<MaskCalculatorAngular, bases<MaskCalculatorBase>>("MaskCalculatorAngular",
                                                             no_init)
      .def(init<Experiment, mat3<double>, double, mat3<double>>());
  }

}}}  // namespace dials::algorithms::boost_python
