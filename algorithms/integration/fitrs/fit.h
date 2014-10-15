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

#ifndef DIALS_ALGORITHMS_INTEGRATION_FITRS_FIT_H
#define DIALS_ALGORITHMS_INTEGRATION_FITRS_FIT_H

#include <iostream>
#include <iomanip>
#include <scitbx/array_family/tiny_types.h>
#include <dxtbx/model/beam.h>
#include <dxtbx/model/detector.h>
#include <dxtbx/model/goniometer.h>
#include <dxtbx/model/scan.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/array_family/reflection_table.h>
#include <dials/algorithms/integration/profile/fitting.h>
#include <dials/algorithms/reflection_basis/coordinate_system.h>
#include <dials/algorithms/reflection_basis/transform.h>
#include <dials/algorithms/polygon/clip/clip.h>
#include <dials/algorithms/polygon/spatial_interpolation.h>
#include <dials/algorithms/polygon/area.h>
#include <dials/algorithms/integration/profile/reference_learner.h>
#include <dials/algorithms/integration/profile/grid_sampler.h>


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
  using dials::algorithms::reflection_basis::transform::TransformSpec;
  using dials::algorithms::reflection_basis::transform::Forward;
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
         double sigma_b,
         double sigma_m,
         double n_sigma)
      : beam_(beam),
        detector_(detector),
        goniometer_(goniometer),
        scan_(scan),
        sigma_b_(sigma_b),
        sigma_m_(sigma_m),
        n_sigma_(n_sigma) {
      DIALS_ASSERT(sigma_b_ > 0);
      DIALS_ASSERT(sigma_m_ > 0);
      DIALS_ASSERT(n_sigma_ > 0);
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

    const double sigma_b() const {
      return sigma_b_;
    }

    const double sigma_m() const {
      return sigma_m_;
    }

    const double n_sigma() const {
      return n_sigma_;
    }

  private:

    Beam beam_;
    Detector detector_;
    Goniometer goniometer_;
    Scan scan_;
    double sigma_b_;
    double sigma_m_;
    double n_sigma_;
  };


  class ProfileLearnerIface {
  public:

    typedef double float_type;
    typedef af::const_ref< float_type, af::c_grid<3> > profile_type;
    typedef af::const_ref< bool, af::c_grid<3> > mask_type;

    virtual
    ~ProfileLearnerIface() {}

    virtual
    bool add(const Shoebox<> &sbox, vec3<double> s1, vec3<double> xyz, double phi) = 0;

    virtual
    void finalize() = 0;

    virtual
    profile_type data(std::size_t index) const = 0;

    virtual
    std::size_t size() const = 0;

    virtual
    std::size_t count() const = 0;

    virtual
    profile_type get(vec3<double> xyz) const = 0;

    virtual
    mask_type get_mask(vec3<double> xyz) const = 0;
  };

  class SingleProfileLearner : public ProfileLearnerIface {
  public:

    SingleProfileLearner(
        const Spec &spec,
        std::size_t grid_size,
        double threshold)
      : learner_(
          GridSampler(
            int3(
              spec.detector()[0].get_image_size()[0],
              spec.detector()[0].get_image_size()[1],
              spec.scan().get_num_images()),
            int3(1, 1, 1)),
          int3(
            2 * grid_size + 1,
            2 * grid_size + 1,
            2 * grid_size + 1),
          threshold),
        spec_(
            spec.beam(),
            spec.detector(),
            spec.goniometer(),
            spec.scan(),
            spec.sigma_b(),
            spec.sigma_m(),
            5.0,
            grid_size),
        count_(0),
        finalized_(false) {
      DIALS_ASSERT(spec.detector().size() == 1);
    }

    /**
     * Add a reflection to the reference profile.
     * @param sbox The shoebox
     * @param s1 The diffracted beam vector.
     * @param phi The rotation angle.
     */
    virtual
    bool add(const Shoebox<> &sbox, vec3<double> s1, vec3<double> xyz, double phi) {
      DIALS_ASSERT(!finalized_);
      Forward<double> transform(spec_, s1, phi, sbox);
      learner_.add(transform.profile().const_ref(), xyz);
      count_++;
      return true;
    }

    /**
     * Finialize the profile.
     */
    virtual
    void finalize() {
      DIALS_ASSERT(!finalized_);
      learner_.finalize();
      finalized_ = true;
    }

    /**
     * @returns The profile.
     */
    virtual
    profile_type data(std::size_t index) const {
      DIALS_ASSERT(finalized_);
      return learner_.locate().profile_ref(index);
    }

    virtual
    std::size_t size() const {
      DIALS_ASSERT(finalized_);
      return learner_.locate().size();
    }

    /**
     * @returns The number of reflections used to compute the profile.
     */
    virtual
    std::size_t count() const {
      return count_;
    }

    /**
     * Get the profile for the given reflection.
     * @param The panel number
     * @param s1 The diffracted beam vector
     * @param phi The rotation angle
     * @param bbox The bounding box
     * @returns The profile for the reflection
     */
    virtual
    profile_type get(vec3<double> xyz) const {
      DIALS_ASSERT(finalized_);
      return learner_.locate().profile_ref(xyz);
    }

    virtual
    mask_type get_mask(vec3<double> xyz) const {
      DIALS_ASSERT(finalized_);
      return learner_.locate().mask_ref(xyz);
    }

  private:

    ReferenceLearnerNew<GridSampler> learner_;
    TransformSpec<double> spec_;
    std::size_t count_;
    bool finalized_;
  };

  class GridProfileLearner : public ProfileLearnerIface {
  public:

    GridProfileLearner(
        const Spec &spec,
        std::size_t grid_size,
        double threshold)
      : learner_(
          GridSampler(
            int3(
              spec.detector()[0].get_image_size()[0],
              spec.detector()[0].get_image_size()[1],
              spec.scan().get_num_images()),
            int3(3, 3, 1)),
          int3(
            2 * grid_size + 1,
            2 * grid_size + 1,
            2 * grid_size + 1),
          threshold),
        spec_(
            spec.beam(),
            spec.detector(),
            spec.goniometer(),
            spec.scan(),
            spec.sigma_b(),
            spec.sigma_m(),
            5.0,
            grid_size),
        count_(0),
        finalized_(false) {
      DIALS_ASSERT(spec.detector().size() == 1);
    }

    /**
     * Add a reflection to the reference profile.
     * @param sbox The shoebox
     * @param s1 The diffracted beam vector.
     * @param phi The rotation angle.
     */
    virtual
    bool add(const Shoebox<> &sbox, vec3<double> s1, vec3<double> xyz, double phi) {
      DIALS_ASSERT(!finalized_);
      Forward<double> transform(spec_, s1, phi, sbox);
      learner_.add(transform.profile().const_ref(), xyz);
      count_++;
      return true;
    }

    /**
     * Finialize the profile.
     */
    virtual
    void finalize() {
      DIALS_ASSERT(!finalized_);
      learner_.finalize();
      finalized_ = true;
    }

    /**
     * @returns The profile.
     */
    virtual
    profile_type data(std::size_t index) const {
      DIALS_ASSERT(finalized_);
      return learner_.locate().profile_ref(index);
    }

    virtual
    std::size_t size() const {
      DIALS_ASSERT(finalized_);
      return learner_.locate().size();
    }

    /**
     * @returns The number of reflections used to compute the profile.
     */
    virtual
    std::size_t count() const {
      return count_;
    }

    /**
     * Get the profile for the given reflection.
     * @param The panel number
     * @param s1 The diffracted beam vector
     * @param phi The rotation angle
     * @param bbox The bounding box
     * @returns The profile for the reflection
     */
    virtual
    profile_type get(vec3<double> xyz) const {
      DIALS_ASSERT(finalized_);
      return learner_.locate().profile_ref(xyz);
    }

    virtual
    mask_type get_mask(vec3<double> xyz) const {
      DIALS_ASSERT(finalized_);
      return learner_.locate().mask_ref(xyz);
    }

  private:

    ReferenceLearnerNew<GridSampler> learner_;
    TransformSpec<double> spec_;
    std::size_t count_;
    bool finalized_;
  };

  /**
   * A class to compute reference profiles for each experiment.
   */
  class MultiExpProfileLearner {
  public:

    typedef ProfileLearnerIface learner_type;
    typedef learner_type::float_type float_type;
    typedef learner_type::profile_type profile_type;
    typedef learner_type::mask_type mask_type;

    /**
     * Compute the reference profiles.
     * @param spec The list of experiment specifications
     * @param data The list of reflections.
     */
    MultiExpProfileLearner(
        af::const_ref<Spec> spec,
        af::reflection_table data,
        std::size_t grid_size,
        double threshold,
        bool single_reference) {

      // Check the input
      DIALS_ASSERT(spec.size() > 0);
      DIALS_ASSERT(data.size() > 0);
      DIALS_ASSERT(grid_size > 0);

      // Create the array of learners
      learner_.reserve(spec.size());
      for (std::size_t i = 0; i < spec.size(); ++i) {
        if (single_reference || spec[i].detector().size() > 1) {
          learner_.push_back(
              boost::shared_ptr<learner_type>(
                new SingleProfileLearner(
                  spec[i],
                  grid_size,
                  threshold)));
        } else {
          learner_.push_back(
              boost::shared_ptr<learner_type>(
                new GridProfileLearner(
                  spec[i],
                  grid_size,
                  threshold)));
        }
      }

      // Check the input contains expected fields
      DIALS_ASSERT(data.size() > 0);
      DIALS_ASSERT(data.is_consistent());
      DIALS_ASSERT(data.contains("id"));
      DIALS_ASSERT(data.contains("shoebox"));
      DIALS_ASSERT(data.contains("s1"));
      DIALS_ASSERT(data.contains("xyzcal.mm"));
      DIALS_ASSERT(data.contains("xyzcal.px"));
      DIALS_ASSERT(data.contains("flags"));

      // Get the data we need
      af::const_ref<std::size_t>    id      = data["id"];
      af::const_ref< Shoebox<> >    shoebox = data["shoebox"];
      af::const_ref< vec3<double> > s1      = data["s1"];
      af::const_ref< vec3<double> > xyzpx   = data["xyzcal.px"];
      af::const_ref< vec3<double> > xyzmm   = data["xyzcal.mm"];
      af::ref<std::size_t>          flags   = data["flags"];

      // Loop through all the reflections
      for (std::size_t i = 0; i < data.size(); ++i) {
        DIALS_ASSERT(id[i] < learner_.size());
        if (flags[i] & ReferenceSpot) {
          learner_[id[i]]->add(shoebox[i], s1[i], xyzpx[i], xyzmm[i][2]);
        }
      }

      // Finalize the profiles
      for (std::size_t i = 0; i < learner_.size(); ++i) {
        learner_[i]->finalize();
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
        vec3<double> xyz) const {
      DIALS_ASSERT(id < learner_.size());
      return learner_[id]->get(xyz);
    }

    mask_type get_mask(
        std::size_t id,
        vec3<double> xyz) const {
      DIALS_ASSERT(id < learner_.size());
      return learner_[id]->get_mask(xyz);
    }

    /**
     * Get the reference profile.
     * @param id The experiment ID
     * @returns The reference profile.
     */
    profile_type data(std::size_t id, std::size_t index) {
      DIALS_ASSERT(id < learner_.size());
      return learner_[id]->data(index);
    }

    /**
     * Get the number of profiles used
     * @param id The experiment ID
     * @returns The number of profiles used.
     */
    std::size_t count(std::size_t id) {
      DIALS_ASSERT(id < learner_.size());
      return learner_[id]->count();
    }

    /**
     * Return the number of references.
     */
    std::size_t size() const {
      return learner_.size();
    }

    std::size_t single_size(std::size_t id) const {
      DIALS_ASSERT(id < learner_.size());
      return learner_[id]->size();
    }

  private:

    std::vector<boost::shared_ptr<learner_type> > learner_;
  };


  /**
   * A class to perform reciprocal space profile fitting.
   */
  class ReciprocalSpaceProfileFitting {
  public:

    typedef Shoebox<>::float_type float_type;
    typedef MultiExpProfileLearner reference_learner_type;
    typedef reference_learner_type::profile_type profile_type;
    typedef reference_learner_type::mask_type mask_type;
    typedef ProfileFitting<double> fitting_type;

    ReciprocalSpaceProfileFitting(std::size_t grid_size,
                                  double threshold,
                                  bool single_reference)
        : grid_size_(grid_size),
          threshold_(threshold),
          single_reference_(single_reference) {
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
      DIALS_ASSERT(data.contains("xyzcal.px"));
      DIALS_ASSERT(data.contains("xyzcal.mm"));
      DIALS_ASSERT(data.contains("flags"));

      // Get the data we need
      af::const_ref<std::size_t>    id      = data["id"];
      af::const_ref< Shoebox<> >    shoebox = data["shoebox"];
      af::const_ref< vec3<double> > s1      = data["s1"];
      af::const_ref< vec3<double> > xyzpx  = data["xyzcal.px"];
      af::const_ref< vec3<double> > xyzmm  = data["xyzcal.mm"];
      af::ref<std::size_t>          flags   = data["flags"];

      // Remove any reflections which have invalid pixels or are not fully
      // recorded from integration.
      for (std::size_t i = 0; i < data.size(); ++i) {
        Shoebox<> sbox = shoebox[i];
        std::size_t panel = sbox.panel;
        int6 bbox = sbox.bbox;
        DIALS_ASSERT(id[i] < spec_.size());
        int2 image_size = spec_[id[i]].detector()[panel].get_image_size();
        DIALS_ASSERT(sbox.is_consistent());
        if (bbox[0] < 0 || bbox[1] > image_size[0] ||
            bbox[2] < 0 || bbox[3] > image_size[1]) {
          flags[i] |= DontIntegrate;
          flags[i] &= ~ReferenceSpot;
          continue;
        }
        for (std::size_t j = 0; j < sbox.mask.size(); ++j) {
          if (sbox.mask[j] & Foreground && !(sbox.mask[j] & Valid)) {
            flags[i] |= DontIntegrate;
            flags[i] &= ~ReferenceSpot;
            break;
          }
        }
      }

      // Compute the reference profile
      reference_learner_type reference(
          spec_.const_ref(),
          data,
          grid_size_,
          threshold_,
          single_reference_);

      // Create a list of transform specs
      std::vector< TransformSpec<> > transform_spec;
      for (std::size_t i = 0; i < spec_.size(); ++i) {
        transform_spec.push_back(
            TransformSpec<>(
              spec_[i].beam(),
              spec_[i].detector(),
              spec_[i].goniometer(),
              spec_[i].scan(),
              spec_[i].sigma_b(),
              spec_[i].sigma_m(),
              5.0,
              grid_size_));
      }

      // Get the new columns to set
      af::ref<double> intensity   = data["intensity.prf.value"];
      af::ref<double> variance    = data["intensity.prf.variance"];
      af::ref<double> correlation = data["profile_correlation"];

      // Do the profile fitting for all reflections
      for (std::size_t i = 0; i < data.size(); ++i) {

        // Initialise some stuff
        intensity[i] = 0;
        variance[i] = -1;
        correlation[i] = 0;
        flags[i] &= ~IntegratedPrf;

        // Integrate
        if (!(flags[i] & DontIntegrate)) {

          // Get the shoebox
          const Shoebox<> &sbox = shoebox[i];
          DIALS_ASSERT(sbox.is_consistent());
          DIALS_ASSERT(id[i] < spec_.size());

          // Do the transform
          Forward<> transform(
              transform_spec[id[i]],
              s1[i],
              xyzmm[i][2],
              sbox);

          // Get the profile for a given reflection
          profile_type c = transform.profile().const_ref();
          profile_type b = transform.background().const_ref();
          profile_type p = reference.get(id[i], xyzpx[i]);
          mask_type    m = reference.get_mask(id[i], xyzpx[i]);

          // Perform the profile fit
          try {
            fitting_type fit(p, m, c, b, 1e-3, 10);
            DIALS_ASSERT(fit.niter() < 10);

            // Set the data in the reflection
            intensity[i]   = fit.intensity();
            variance[i]    = fit.variance();
            correlation[i] = fit.correlation();

            // Set the integrated flag
            flags[i] |= IntegratedPrf;
          } catch (dials::error) {
            continue;
          }
        }
      }

      // Return the reference learner
      return reference;
    }

  private:

    af::shared<Spec> spec_;
    std::size_t grid_size_;
    double threshold_;
    bool single_reference_;
  };

}} // namespace dials::algorithms


#endif // DIALS_ALGORITHMS_INTEGRATION_FITRS_FIT_H
