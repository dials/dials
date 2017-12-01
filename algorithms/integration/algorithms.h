/*
 * algorithms.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */

#ifndef DIALS_ALGORITHMS_INTEGRATION_ALGORITHMS_H
#define DIALS_ALGORITHMS_INTEGRATION_ALGORITHMS_H

#include <numeric>

#include <dxtbx/model/beam.h>
#include <dxtbx/model/detector.h>
#include <dxtbx/model/scan.h>
#include <dxtbx/model/crystal.h>
#include <dials/algorithms/profile_model/gaussian_rs/mask_calculator.h>

#include <dials/algorithms/background/simple/creator.h>
#include <dials/algorithms/background/glm/creator.h>
#include <dials/algorithms/background/gmodel/creator.h>

#include <dials/algorithms/profile_model/modeller/circle_sampler.h>
#include <dials/algorithms/profile_model/gaussian_rs/transform/transform.h>
#include <dials/algorithms/integration/fit/fitting.h>
#include <dials/algorithms/integration/interfaces.h>
#include <dials/algorithms/spot_prediction/pixel_to_miller_index.h>
#include <dials/model/data/mask_code.h>
#include <dials/model/data/shoebox.h>

namespace dials { namespace algorithms {

  using dxtbx::model::Beam;
  using dxtbx::model::Detector;
  using dxtbx::model::Goniometer;
  using dxtbx::model::Scan;

  using dials::algorithms::profile_model::gaussian_rs::MaskCalculator3D;
  using dials::algorithms::profile_model::gaussian_rs::transform::TransformSpec;
  using dials::algorithms::profile_model::gaussian_rs::CoordinateSystem;
  using dials::algorithms::profile_model::gaussian_rs::transform::TransformReverseNoModel;
  using dials::algorithms::profile_model::gaussian_rs::transform::TransformReverse;
  using dials::algorithms::profile_model::gaussian_rs::transform::TransformForward;
  using dials::algorithms::background::SimpleBackgroundCreator;
  using dials::algorithms::GLMBackgroundCreator;
  using dials::algorithms::GModelBackgroundCreator;


  /**
   * A class to calculate the reflection mask
   */
  class GaussianRSMaskCalculator : public MaskCalculatorIface {
  public:

    /**
     * Init the algorithm
     * @param beam The beam model
     * @param detector The detector model
     * @param gonio The goniometer model
     * @param scan The scan model
     * @param delta_b The beam divergence
     * @param delta_m The mosaicity
     */
    GaussianRSMaskCalculator(const BeamBase &beam,
                             const Detector &detector,
                             const Goniometer &gonio,
                             const Scan &scan,
                             double delta_b,
                             double delta_m)
      : func_(beam,
              detector,
              gonio,
              scan,
              delta_b,
              delta_m) {}

    /**
     * Compute the mask for a single reflection
     * @param reflection The reflection object
     * @param adjacent Is this an adjacent relfection?
     */
    virtual void operator()(af::Reflection &reflection, bool adjacent=false) const {
      func_.single(
        reflection.get< Shoebox<> >("shoebox"),
        reflection.get< vec3<double> >("s1"),
        reflection.get< vec3<double> >("xyzcal.px")[2],
        reflection.get< std::size_t >("panel"),
        adjacent);
    }

  protected:

    MaskCalculator3D func_;

  };


  /**
   * A class to calculate the reflection mask
   */
  class GaussianRSMultiCrystalMaskCalculator : public MaskCalculatorIface {
  public:

    /**
     * Add a mask calculator to the list
     * @param alg The mask calculator
     */
    void append(const GaussianRSMaskCalculator &alg) {
      algorithms_.push_back(alg);
    }

    /**
     * Compute the mask for a single reflection
     * @param reflection The reflection object
     * @param adjacent Is this an adjacent relfection?
     */
    virtual void operator()(af::Reflection &reflection, bool adjacent=false) const {
      int index = reflection.get<int>("id");
      DIALS_ASSERT(index >= 0 && index < algorithms_.size());
      algorithms_[index](reflection, adjacent);
    }

  protected:

    std::vector<GaussianRSMaskCalculator> algorithms_;

  };


  /**
   * A class to do simple background determination
   */
  class SimpleBackgroundCalculator : public BackgroundCalculatorIface {
  public:

    /**
     * Init the algorithm
     * @param modeller The modeller algorithm
     * @param rejector The outlier rejector
     * @param min_pixels The min number of pixels
     */
    SimpleBackgroundCalculator(
          boost::shared_ptr<background::Modeller> modeller,
          boost::shared_ptr<background::OutlierRejector> rejector,
          std::size_t min_pixels)
      : creator_(
          modeller,
          rejector,
          min_pixels) {}


    /**
     * Compute the background
     * @param reflection The reflection object
     */
    virtual void operator()(af::Reflection &reflection) const {
      creator_(reflection.get< Shoebox<> >("shoebox"));
    }

  protected:

    SimpleBackgroundCreator creator_;
  };


  /**
   * A class to do GLM background determination
   */
  class GLMBackgroundCalculator : public BackgroundCalculatorIface {
  public:

    /**
     * Init the algorithms
     * @param model The model to use
     * @param tuning_constant The robust tuning constant
     * @param max_iter The maximum number of iterations
     * @param min_pixels The minimum number of pixels
     */
    GLMBackgroundCalculator(
          GLMBackgroundCreator::Model model,
          double tuning_constant,
          std::size_t max_iter,
          std::size_t min_pixels)
      : creator_(
          model,
          tuning_constant,
          max_iter,
          min_pixels) {}

    /**
     * Compute the background
     * @param reflection The reflection object
     */
    virtual void operator()(af::Reflection &reflection) const {
      creator_.single(reflection.get< Shoebox<> >("shoebox"));
    }

  protected:

    GLMBackgroundCreator creator_;
  };


  /**
   * A class to do global model background determination
   */
  class GModelBackgroundCalculator : public BackgroundCalculatorIface {
  public:

    /**
     * Init the algorithms
     * @param model The model to use
     * @param robust Do robust estimation
     * @param tuning_constant The robust tuning constant
     * @param max_iter The maximum number of iterations
     * @param min_pixels The minimum number of pixels
     */
    GModelBackgroundCalculator(
          boost::shared_ptr<BackgroundModel> model,
          bool robust,
          double tuning_constant,
          std::size_t max_iter,
          std::size_t min_pixels)
      : creator_(
          model,
          robust,
          tuning_constant,
          max_iter,
          min_pixels) {}

    /**
     * Compute the background
     * @param reflection The reflection object
     */
    virtual void operator()(af::Reflection &reflection) const {
      creator_.single(reflection.get< Shoebox<> >("shoebox"));
    }

  protected:

    GModelBackgroundCreator creator_;
  };


  /**
   * A class implementing the intentity calculator interface but doing nothing
   */
  class NullIntensityCalculator : public IntensityCalculatorIface {
  public:

    /**
     * Do nothing
     */
    virtual void operator()(af::Reflection &reflection,
                            const std::vector<af::Reflection> &adjacent_reflections) const {

    }

  };


  /**
   * A class to hold the reference data
   */
  class ReferenceProfileData {
  public:

    typedef af::versa< double, af::c_grid<3> > data_type;
    typedef af::versa< bool, af::c_grid<3> > mask_type;

    /**
     * Append the reference data
     * @param data The data
     * @param mask The mask
     */
    void append(af::const_ref< double, af::c_grid<3> > data,
                af::const_ref< bool, af::c_grid<3> > mask) {
      DIALS_ASSERT(data.accessor().all_eq(mask.accessor()));
      data_type d(data.accessor());
      mask_type m(mask.accessor());
      std::copy(data.begin(), data.end(), d.begin());
      std::copy(mask.begin(), mask.end(), m.begin());
      data_.push_back(d);
      mask_.push_back(m);
    }

    /**
     * Get the data at the index
     * @param index The reference profile index
     * @returns The reference profile
     */
    af::const_ref< double, af::c_grid<3> > data(std::size_t index) const {
      DIALS_ASSERT(index < data_.size());
      return data_[index].const_ref();
    }

    /**
     * Get the mask at the index
     * @param index The reference profile index
     * @returns The reference profile mask
     */
    af::const_ref< bool, af::c_grid<3> > mask(std::size_t index) const {
      DIALS_ASSERT(index < mask_.size());
      return mask_[index].const_ref();
    }

    /**
     * @returns number of reference profiles
     */
    std::size_t size() const {
      DIALS_ASSERT(data_.size() == mask_.size());
      return data_.size();
    }

  protected:

    af::shared<data_type> data_;
    af::shared<mask_type> mask_;
  };


  /**
   * A class to hold all the reference data
   */
  class GaussianRSReferenceProfileData {
  public:

    /**
     * Initialise the algorithm
     * @param reference The reference profiles
     * @param sampler The sampler object
     * @param spec The transform spec
     */
    GaussianRSReferenceProfileData(
        const ReferenceProfileData &reference,
        boost::shared_ptr<SamplerIface> sampler,
        const TransformSpec &spec)
      : reference_(reference),
        sampler_(sampler),
        spec_(spec) {}

    /**
     * Get the reference
     */
    const ReferenceProfileData& reference() const {
      return reference_;
    }

    /**
     * Get the sampler
     */
    boost::shared_ptr<SamplerIface> sampler() const {
      return sampler_;
    }

    /**
     * Get the transform spec
     */
    const TransformSpec& spec() const {
      return spec_;
    }

  protected:

    ReferenceProfileData reference_;
    boost::shared_ptr<SamplerIface> sampler_;
    TransformSpec spec_;

  };


  /**
   * A class to hold data for multiple crystals
   */
  class GaussianRSMultiCrystalReferenceProfileData {
  public:

    /**
     * Add a data spec to the list
     * @param alg The mask calculator
     */
    void append(const GaussianRSReferenceProfileData &spec) {
      spec_list_.push_back(spec);
    }

    /**
     * Get the data spec for the index
     * @param index The index of the experiment
     * @returns The data spec for the experiments
     */
    const GaussianRSReferenceProfileData& operator[](std::size_t index) const {
      DIALS_ASSERT(index < spec_list_.size());
      return spec_list_[index];
    }

    /**
     * @returns number of data structs
     */
    std::size_t size() const {
      return spec_list_.size();
    }

  protected:

    std::vector<GaussianRSReferenceProfileData> spec_list_;

  };


  /**
   * Interface to different GaussianRS intensity algorithms
   */
  class GaussianRSIntensityCalculatorAlgorithm {
  public:

    virtual ~GaussianRSIntensityCalculatorAlgorithm() {}

    virtual
    void exec(
          af::Reflection &reflection,
          const std::vector<af::Reflection> &adjacent_reflections) const = 0;

  };


  /**
   * A class to implement standard reciprocal space profile fitting
   */
  class GaussianRSReciprocalSpaceIntensityCalculator
    : public GaussianRSIntensityCalculatorAlgorithm {
  public:

    /**
     * Initialise the algorithm
     * @param data The data to do profile fitting
     */
    GaussianRSReciprocalSpaceIntensityCalculator(
        const GaussianRSMultiCrystalReferenceProfileData &data)
      : data_spec_(data) {}

    /**
     * Compute the intensity
     * @param reflection The reflection object
     * @param adjacent_reflections The adjacent reflections
     */
    virtual
    void exec(
          af::Reflection &reflection,
          const std::vector<af::Reflection> &adjacent_reflections) const {

      // Typedefs
      typedef af::const_ref< double, af::c_grid<3> > data_const_reference;
      typedef af::const_ref< bool, af::c_grid<3> >   mask_const_reference;

      // Get the index of the reflection
      int experiment_id = reflection.get< int >("id");
      DIALS_ASSERT(experiment_id >= 0);
      const GaussianRSReferenceProfileData &data_spec = data_spec_[experiment_id];

      // Get the reflection flags and Unset profile fitting
      reflection["flags"] = reflection.get<std::size_t>("flags") & ~af::IntegratedPrf;

      // Get some stuff from the reflection
      vec3<double> s1 = reflection.get< vec3<double> >("s1");
      double phi = reflection.get< vec3<double> >("xyzcal.mm")[2];
      vec3<double> xyz = reflection.get< vec3<double> >("xyzcal.px");
      Shoebox<> sbox = reflection.get< Shoebox<> >("shoebox");
      DIALS_ASSERT(sbox.is_consistent());

      // Check the shoebox
      if (check(sbox, data_spec.spec().detector()) == false) {
        throw DIALS_ERROR("Invalid pixels in Shoebox");
      }

      // Create the coordinate system
      vec3<double> m2 = data_spec.spec().goniometer().get_rotation_axis();
      vec3<double> s0 = data_spec.spec().beam()->get_s0();
      CoordinateSystem cs(m2, s0, s1, phi);

      // Get the reference profiles
      std::size_t index = data_spec.sampler()->nearest(sbox.panel, xyz);
      data_const_reference reference_data = data_spec.reference().data(index);
      mask_const_reference reference_mask = data_spec.reference().mask(index);

      // Create the data array
      af::versa< double, af::c_grid<3> > data(sbox.data.accessor());
      std::copy(
          sbox.data.begin(),
          sbox.data.end(),
          data.begin());

      // Create the background array
      af::versa< double, af::c_grid<3> > background(sbox.background.accessor());
      std::copy(
          sbox.background.begin(),
          sbox.background.end(),
          background.begin());

      // Create the mask array
      af::versa< bool, af::c_grid<3> > mask(sbox.mask.accessor());
      copy_mask(
          sbox.mask.begin(),
          sbox.mask.end(),
          mask.begin());

      // Compute the transform
      TransformForward<double> transform(
          data_spec.spec(),
          cs,
          sbox.bbox,
          sbox.panel,
          data.const_ref(),
          background.const_ref(),
          mask.const_ref());

      // Get the transformed shoebox
      data_const_reference transformed_data = transform.profile().const_ref();
      data_const_reference transformed_bgrd = transform.background().const_ref();
      mask_const_reference transformed_mask = transform.mask().const_ref();
      af::versa< bool, af::c_grid<3> > final_mask(transformed_mask.accessor());
      DIALS_ASSERT(reference_mask.size() == transformed_mask.size());
      std::size_t mask_count = 0;
      for (std::size_t j = 0; j < final_mask.size(); ++j) {
        final_mask[j] = reference_mask[j] && transformed_mask[j];
        if (final_mask[j]) {
          mask_count++;
        }
      }
      if (mask_count == 0) {
        throw DIALS_ERROR("No pixels mapped to reciprocal space grid");
      }

      // Do the profile fitting
      ProfileFitter<double> fit(
          transformed_data,
          transformed_bgrd,
          final_mask.const_ref(),
          reference_data,
          1e-3,
          100);
      DIALS_ASSERT(fit.niter() < 100);

      // Set the integrated data and the integrated flag
      reflection["intensity.prf.value"] = fit.intensity()[0];
      reflection["intensity.prf.variance"] = fit.variance()[0];
      reflection["intensity.prf.correlation"] = fit.correlation();
      reflection["flags"] = reflection.get<std::size_t>("flags") | af::IntegratedPrf;
    }

  protected:

    /**
     * Copy the mask while checking the mask codes
     */
    template <typename InputIterator, typename OutputIterator>
    void copy_mask(InputIterator first,
                   InputIterator last,
                   OutputIterator out) const {
      std::size_t mask_code = Valid | Foreground;
      for (InputIterator it = first; it != last; ++it) {
        *out++ = ((*it & mask_code) == mask_code && (*it & Overlapped) == 0);
      }
    }

    /**
     * @returns True/False if the shoebox is valid
     */
    bool check(const Shoebox<> &sbox, const Detector &detector) const {

      // Check if the bounding box is in the image
      bool bbox_valid =
        sbox.bbox[0] >= 0 &&
        sbox.bbox[2] >= 0 &&
        sbox.bbox[1] <= detector[sbox.panel].get_image_size()[0] &&
        sbox.bbox[3] <= detector[sbox.panel].get_image_size()[1];

      // Check if all pixels are valid
      bool pixels_valid = true;
      for (std::size_t i = 0; i < sbox.mask.size(); ++i) {
        if ((sbox.mask[i] & Foreground) && (
            ((sbox.mask[i] & Valid) == 0) ||
            (sbox.mask[i] & Overlapped))) {
          pixels_valid = false;
          break;
        }
      }

      // Return whether to use or not
      return bbox_valid && pixels_valid;
    }

    GaussianRSMultiCrystalReferenceProfileData data_spec_;

  };


  /**
   * Class to implement profile fitting in detector space
   */
  class GaussianRSDetectorSpaceIntensityCalculator
    : public GaussianRSIntensityCalculatorAlgorithm {
  public:

    /**
     * Initialise the algorithm
     * @param data The data to do profile fitting
     */
    GaussianRSDetectorSpaceIntensityCalculator(
        const GaussianRSMultiCrystalReferenceProfileData &data)
      : data_spec_(data) {}

    /**
     * Compute the intensity
     * @param reflection The reflection object
     * @param adjacent_reflections The adjacent reflections
     */
    virtual
    void exec(
          af::Reflection &reflection,
          const std::vector<af::Reflection> &adjacent_reflections) const {

      // Typedefs
      typedef af::const_ref< double, af::c_grid<3> > data_const_reference;
      typedef af::const_ref< bool, af::c_grid<3> > mask_const_reference;

      // Get the index of the reflection
      int experiment_id = reflection.get< int >("id");
      DIALS_ASSERT(experiment_id >= 0);
      const GaussianRSReferenceProfileData &data_spec = data_spec_[experiment_id];

      // Get the reflection flags and Unset profile fitting
      reflection["flags"] = reflection.get<std::size_t>("flags") & ~af::IntegratedPrf;

      // Get some data
      vec3<double> s1 = reflection.get< vec3<double> >("s1");
      double phi = reflection.get< vec3<double> >("xyzcal.mm")[2];
      vec3<double> xyz = reflection.get< vec3<double> >("xyzcal.px");
      Shoebox<> sbox = reflection.get< Shoebox<> >("shoebox");
      DIALS_ASSERT(sbox.is_consistent());

      // Get the reference profiles
      std::size_t index = data_spec.sampler()->nearest(sbox.panel, xyz);
      data_const_reference transformed_reference_data = data_spec.reference().data(index);
      //mask_const_reference transformed_reference_mask = data_spec.reference().mask(index);

      // Create the coordinate system
      vec3<double> m2 = data_spec.spec().goniometer().get_rotation_axis();
      vec3<double> s0 = data_spec.spec().beam()->get_s0();
      CoordinateSystem cs(m2, s0, s1, phi);

      // Compute the transform
      TransformReverse transform(
          data_spec.spec(),
          cs,
          sbox.bbox,
          sbox.panel,
          transformed_reference_data);

      // Get the transformed shoebox
      data_const_reference reference_data = transform.profile().const_ref();

      // Create the data array
      af::versa< double, af::c_grid<3> > data(sbox.data.accessor());
      std::copy(
          sbox.data.begin(),
          sbox.data.end(),
          data.begin());

      // Create the background array
      af::versa< double, af::c_grid<3> > bgrd(sbox.background.accessor());
      std::copy(
          sbox.background.begin(),
          sbox.background.end(),
          bgrd.begin());

      // Create the mask array
      af::versa< bool, af::c_grid<3> > mask(sbox.mask.accessor());
      copy_mask(
        sbox.mask.begin(),
        sbox.mask.end(),
        mask.begin());

      // Compute partiality
      double partiality = compute_partiality(reference_data, mask.const_ref());
      DIALS_ASSERT(partiality >= 0 && partiality <= 1.0);

      // Exit if partiality too small
      if (partiality < 0.01) {
        throw DIALS_ERROR("Partiality too small to fit");
      }

      // Do the profile fitting
      ProfileFitter<double> fit(
          data.const_ref(),
          bgrd.const_ref(),
          mask.const_ref(),
          reference_data,
          1e-3, 100);
      DIALS_ASSERT(fit.niter() < 100);

      // Set the integrated values and flag
      reflection["intensity.prf.value"] = fit.intensity()[0];
      reflection["intensity.prf.variance"] = fit.variance()[0];
      reflection["intensity.prf.correlation"] = fit.correlation();
      reflection["partiality_old"] = reflection.get<double>("partiality");
      reflection["partiality"] = partiality;
      reflection["flags"] = reflection.get<std::size_t>("flags") | af::IntegratedPrf;
    }

  protected:

    /**
     * Compute partiality
     */
    double compute_partiality(
        const af::const_ref< double, af::c_grid<3> > &data,
        const af::const_ref< bool, af::c_grid<3> > &mask) const {
      DIALS_ASSERT(data.size() == mask.size());
      double partiality = 0.0;
      for (std::size_t i = 0; i < data.size(); ++i) {
        if (mask[i]) {
          DIALS_ASSERT(data[i] >= 0);
          partiality += data[i];
        }
      }
      return partiality;
    }

    /**
     * Copy the mask while checking the mask codes
     */
    template <typename InputIterator, typename OutputIterator>
    void copy_mask(InputIterator first,
                   InputIterator last,
                   OutputIterator out) const {
      std::size_t mask_code = Valid | Foreground;
      for (InputIterator it = first; it != last; ++it) {
        *out++ = ((*it & mask_code) == mask_code && (*it & Overlapped) == 0);
      }
    }

    GaussianRSMultiCrystalReferenceProfileData data_spec_;

  };


  /**
   * Class to implement profile fitting with deconvolution
   */
  class GaussianRSDetectorSpaceWithDeconvolutionIntensityCalculator
    : public GaussianRSDetectorSpaceIntensityCalculator {
  public:

    /**
     * Initialise the algorithm
     * @param data The data to do profile fitting
     */
    GaussianRSDetectorSpaceWithDeconvolutionIntensityCalculator(
        const GaussianRSMultiCrystalReferenceProfileData &data)
      : GaussianRSDetectorSpaceIntensityCalculator(data) {}

    /**
     * Compute the intensity
     * @param reflection The reflection object
     * @param adjacent_reflections The adjacent reflections
     */
    virtual
    void exec(
          af::Reflection &reflection,
          const std::vector<af::Reflection> &adjacent_reflections) const {

      // Typedefs
      typedef af::const_ref< double, af::c_grid<3> > data_const_reference;
      typedef af::const_ref< bool, af::c_grid<3> > mask_const_reference;

      // Get the index of the reflection
      int experiment_id = reflection.get< int >("id");
      DIALS_ASSERT(experiment_id >= 0);
      const GaussianRSReferenceProfileData &data_spec = data_spec_[experiment_id];

      // Get the reflection flags and Unset profile fitting
      reflection["flags"] = reflection.get<std::size_t>("flags") & ~af::IntegratedPrf;

      // Get some data
      vec3<double> s1 = reflection.get< vec3<double> >("s1");
      double phi = reflection.get< vec3<double> >("xyzcal.mm")[2];
      vec3<double> xyz = reflection.get< vec3<double> >("xyzcal.px");
      Shoebox<> sbox = reflection.get< Shoebox<> >("shoebox");
      DIALS_ASSERT(sbox.is_consistent());

      // If we have no overlaps then do normal detector space integration
      if (adjacent_reflections.size() == 0 || !is_overlapping(sbox.mask.const_ref())) {
        GaussianRSDetectorSpaceIntensityCalculator::exec(reflection, adjacent_reflections);
        return;
      }

      // Create the data array
      af::versa< double, af::c_grid<3> > data(sbox.data.accessor());
      std::copy(
          sbox.data.begin(),
          sbox.data.end(),
          data.begin());

      // Create the background array
      af::versa< double, af::c_grid<3> > bgrd(sbox.background.accessor());
      std::copy(
          sbox.background.begin(),
          sbox.background.end(),
          bgrd.begin());

      // Create the mask array
      af::versa< bool, af::c_grid<3> > mask(sbox.mask.accessor());
      copy_mask(
        sbox.mask.begin(),
        sbox.mask.end(),
        mask.begin());

      // Get the reference profiles
      std::size_t index = data_spec.sampler()->nearest(sbox.panel, xyz);
      data_const_reference transformed_reference_data = data_spec.reference().data(index);
      //mask_const_reference transformed_reference_mask = data_spec.reference().mask(index);

      // The profile grid with all adjacent reflections
      af::versa< double, af::c_grid<4> > profile(af::c_grid<4>(
            af::tiny<std::size_t,4>(
            adjacent_reflections.size()+1,
            data.accessor()[0],
            data.accessor()[1],
            data.accessor()[2])));

      // Create the coordinate system
      vec3<double> m2 = data_spec.spec().goniometer().get_rotation_axis();
      vec3<double> s0 = data_spec.spec().beam()->get_s0();
      CoordinateSystem cs(m2, s0, s1, phi);

      // Compute the transform
      TransformReverse transform(
          data_spec.spec(),
          cs,
          sbox.bbox,
          sbox.panel,
          transformed_reference_data);

      // Get the transformed shoebox
      data_const_reference reference_data = transform.profile().const_ref();
      std::copy(
          reference_data.begin(),
          reference_data.end(),
          profile.begin());

      // Compute partiality
      double partiality = compute_partiality(reference_data, mask.const_ref());
      DIALS_ASSERT(partiality >= 0 && partiality <= 1.0);

      // Exit if partiality too small
      if (partiality < 0.01) {
        throw DIALS_ERROR("Partiality too small to fit");
      }

      // Add profiles for adjacent reflections
      for (std::size_t j = 0; j < adjacent_reflections.size(); ++j) {

        // Get the index of the reflection
        int experiment_id = adjacent_reflections[j].get< int >("id");
        DIALS_ASSERT(experiment_id >= 0);
        const GaussianRSReferenceProfileData &data_spec2 = data_spec_[experiment_id];

        // Compute coordinate system
        vec3<double> s12 = adjacent_reflections[j].get< vec3<double> >("s1");
        double phi2 = adjacent_reflections[j].get< vec3<double> >("xyzcal.mm")[2];
        vec3<double> m22 = data_spec2.spec().goniometer().get_rotation_axis();
        vec3<double> s02 = data_spec2.spec().beam()->get_s0();
        CoordinateSystem cs2(m22, s02, s12, phi2);

        // Compute the transform
        TransformReverse transform(
            data_spec2.spec(),
            cs2,
            sbox.bbox,
            sbox.panel,
            transformed_reference_data);

        // Get the transformed shoebox
        data_const_reference reference_data2 = transform.profile().const_ref();

        // Copy into the large profile grid
        std::copy(
            reference_data2.begin(),
            reference_data2.end(),
            profile.begin() + (j+1)* reference_data2.size());
      }

      try {

        // Do the profile fitting
        ProfileFitter<double> fit(
            data.const_ref(),
            bgrd.const_ref(),
            mask.const_ref(),
            profile.const_ref(),
            1e-3, 100);
        DIALS_ASSERT(fit.niter() < 100);

        // Set the integrated data and flag
        reflection["intensity.prf.value"] = fit.intensity()[0];
        reflection["intensity.prf.variance"] = fit.variance()[0];
        reflection["intensity.prf.correlation"] = fit.correlation();
        reflection["partiality_old"] = reflection.get<double>("partiality");
        reflection["partiality"] = partiality;
        reflection["flags"] = reflection.get<std::size_t>("flags") | af::IntegratedPrf;

      } catch (std::runtime_error) {

        // This is only thrown if the matrix is singular, in which case try and
        // do profile fitting with no deconvolution
        GaussianRSDetectorSpaceIntensityCalculator::exec(reflection, adjacent_reflections);
      }
    }

  protected:

    /**
     * Check if the reflection foreground is overlapping
     * @param mask The mask array
     * @returns True/False if the shoebox is overlapping
     */
    bool is_overlapping(const af::const_ref< int, af::c_grid<3> > &mask) const {
      bool result = false;
      int mask_code = Valid | Foreground | Overlapped;
      for (std::size_t i =0 ; i < mask.size(); ++i) {
        if ((mask[i] & mask_code) == mask_code) {
          result = true;
          break;
        }
      }
      return result;
    }

    /**
     * Copy the mask while checking the mask codes
     */
    template <typename InputIterator, typename OutputIterator>
    void copy_mask(InputIterator first,
                   InputIterator last,
                   OutputIterator out) const {
      std::size_t mask_code1 = Valid | Foreground;
      std::size_t mask_code2 = Valid | Overlapped;
      for (InputIterator it = first; it != last; ++it) {
        *out++ = ((*it & mask_code1) == mask_code1) || ((*it & mask_code2) == mask_code2);
      }
    }

  };


  /**
   * A class implementing profile fitting algorithms:
   * 1. In reciprocal space
   * 2. In detector space
   * 3. In detector space with deconvolution
   */
  class GaussianRSIntensityCalculator : public IntensityCalculatorIface {
  public:

    /**
     * Initialise the algorithm
     * @param data The reference profiles
     * @param detector_space Perform in detector space
     * @param deconvolution Do spot deconvolution
     */
    GaussianRSIntensityCalculator(
        const GaussianRSMultiCrystalReferenceProfileData &data,
        bool detector_space,
        bool deconvolution) {

      typedef GaussianRSDetectorSpaceWithDeconvolutionIntensityCalculator DSDCAlgorithm;
      typedef GaussianRSDetectorSpaceIntensityCalculator DSAlgorithm;
      typedef GaussianRSReciprocalSpaceIntensityCalculator RSAlgorithm;

      if (deconvolution) {
        DIALS_ASSERT(detector_space);
        algorithm_ = boost::make_shared<DSDCAlgorithm>(data);
      } else if (detector_space) {
        algorithm_ = boost::make_shared<DSAlgorithm>(data);
      } else {
        algorithm_ = boost::make_shared<RSAlgorithm>(data);
      }
    }

    /**
     * Perform the integration
     * @param reflection The reflection object
     * @param adjacent_reflections The adjacent reflections list
     */
    virtual void operator()(af::Reflection &reflection,
                            const std::vector<af::Reflection> &adjacent_reflections) const {
      algorithm_->exec(reflection, adjacent_reflections);
    }

  protected:

    boost::shared_ptr<GaussianRSIntensityCalculatorAlgorithm> algorithm_;

  };

}}

#endif // DIALS_ALGORITHMS_INTEGRATION_ALGORITHMS_H
