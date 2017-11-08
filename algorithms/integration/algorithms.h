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


  class Reference {
  public:

    typedef af::versa< double, af::c_grid<3> > data_type;
    typedef af::versa< bool, af::c_grid<3> > mask_type;

    void append(af::const_ref< double, af::c_grid<3> > data,
                af::const_ref< bool, af::c_grid<3> > mask) {
      data_type d(data.accessor());
      mask_type m(mask.accessor());
      std::copy(data.begin(), data.end(), d.begin());
      std::copy(mask.begin(), mask.end(), m.begin());
      data_.push_back(d);
      mask_.push_back(m);
    }

    af::const_ref< double, af::c_grid<3> > data(std::size_t index) const {
      DIALS_ASSERT(index < data_.size());
      return data_[index].const_ref();
    }

    af::const_ref< bool, af::c_grid<3> > mask(std::size_t index) const {
      DIALS_ASSERT(index < mask_.size());
      return mask_[index].const_ref();
    }

  protected:

    af::shared<data_type> data_;
    af::shared<mask_type> mask_;
  };

  namespace detail {

    struct check_mask_code {
      int mask_code;
      check_mask_code(int code) : mask_code(code) {}
      bool operator()(int a) const {
        return ((a & mask_code) == mask_code && (a & Overlapped) == 0);
      }
    };

    struct check_mask_code2 {
      int mask_code;
      check_mask_code2(int code) : mask_code(code) {}
      bool operator()(int a) const {
        return ((a & mask_code) == mask_code || (a & (Valid | Overlapped)) == (Valid | Overlapped));
      }
    };

    struct check_either_mask_code {
      int mask_code1;
      int mask_code2;
      check_either_mask_code(int code1, int code2) : mask_code1(code1), mask_code2(code2) {}
      bool operator()(int a) const {
        return ((a & mask_code1) == mask_code1) || ((a & mask_code2) == mask_code2);
      }
    };

  }

  class NullIntensityCalculator : public IntensityCalculatorIface {
  public:

    virtual void operator()(af::Reflection &reflection,
                            const std::vector<af::Reflection> &adjacent_reflections) const {

    }

  };

  class IntensityCalculator : public IntensityCalculatorIface {
  public:

    IntensityCalculator(
        const Reference &reference,
        const CircleSampler &sampler,
        const TransformSpec &spec,
        bool detector_space,
        bool deconvolution)
      : reference_(reference),
        sampler_(sampler),
        spec_(spec),
        detector_space_(detector_space),
        deconvolution_(deconvolution) {
      if (deconvolution) {
        DIALS_ASSERT(detector_space);
      }
    }

    virtual void operator()(af::Reflection &reflection,
                            const std::vector<af::Reflection> &adjacent_reflections) const {
      if (detector_space_ == true) {
        if (deconvolution_ == true) {
          integrate_detector_space_with_deconvolution(
              reflection,
              adjacent_reflections);
        } else {
          integrate_detector_space(
              reflection,
              adjacent_reflections);
        }
      } else {
        integrate_reciprocal_space(
            reflection,
            adjacent_reflections);
      }
    }

    void integrate_reciprocal_space(
          af::Reflection &reflection,
          const std::vector<af::Reflection> &adjacent_reflections) const {

      typedef af::const_ref< double, af::c_grid<3> > data_const_reference;
      typedef af::const_ref< bool, af::c_grid<3> > mask_const_reference;

      reflection["flags"] = reflection.get<std::size_t>("flags") & ~af::IntegratedPrf;

      vec3<double> s1 = reflection.get< vec3<double> >("s1");
      double phi = reflection.get< vec3<double> >("xyzcal.mm")[2];
      vec3<double> xyz = reflection.get< vec3<double> >("xyzcal.px");
      Shoebox<> sbox = reflection.get< Shoebox<> >("shoebox");
      DIALS_ASSERT(sbox.is_consistent());

      if (check(sbox) == false) {
        return;
      }

      // Create the coordinate system
      vec3<double> m2 = spec_.goniometer().get_rotation_axis();
      vec3<double> s0 = spec_.beam()->get_s0();
      CoordinateSystem cs(m2, s0, s1, phi);

      // Get the reference profiles
      std::size_t index = sampler_.nearest(sbox.panel, xyz);
      data_const_reference p = reference_.data(index);
      mask_const_reference mask1 = reference_.mask(index);

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
      std::transform(
          sbox.mask.begin(),
          sbox.mask.end(),
          mask.begin(),
          detail::check_mask_code(Valid | Foreground));

      // Compute the transform
      TransformForward<double> transform(
          spec_,
          cs,
          sbox.bbox,
          sbox.panel,
          data.const_ref(),
          background.const_ref(),
          mask.const_ref());

      // Get the transformed shoebox
      data_const_reference c = transform.profile().const_ref();
      data_const_reference b = transform.background().const_ref();
      mask_const_reference mask2 = transform.mask().const_ref();
      af::versa< bool, af::c_grid<3> > m(mask2.accessor());
      DIALS_ASSERT(mask1.size() == mask2.size());
      for (std::size_t j = 0; j < m.size(); ++j) {
        m[j] = mask1[j] && mask2[j];
      }

      ProfileFitter<double> fit(
          c,
          b,
          m.const_ref(),
          p,
          1e-3,
          100);
      DIALS_ASSERT(fit.niter() < 100);

      // Set the integrated flag
      reflection["intensity.prf.value"] = fit.intensity()[0];
      reflection["intensity.prf.variance"] = fit.variance()[0];
      reflection["intensity.prf.correlation"] = fit.correlation();
      reflection["flags"] = reflection.get<std::size_t>("flags") | af::IntegratedPrf;
    }

    void integrate_detector_space(
          af::Reflection &reflection,
          const std::vector<af::Reflection> &adjacent_reflections) const {

      typedef af::const_ref< double, af::c_grid<3> > data_const_reference;
      typedef af::const_ref< bool, af::c_grid<3> > mask_const_reference;

      reflection["flags"] = reflection.get<std::size_t>("flags") & ~af::IntegratedPrf;

      // Get some data
      vec3<double> s1 = reflection.get< vec3<double> >("s1");
      double phi = reflection.get< vec3<double> >("xyzcal.mm")[2];
      vec3<double> xyz = reflection.get< vec3<double> >("xyzcal.px");
      Shoebox<> sbox = reflection.get< Shoebox<> >("shoebox");
      std::size_t flags = reflection.get<std::size_t>("flags");

      DIALS_ASSERT(sbox.is_consistent());


      // Check if we want to use this reflection
      if (check2(flags, sbox) == false) {
        return;
      }

      // Get the reference profiles
      std::size_t index = sampler_.nearest(sbox.panel, xyz);
      data_const_reference d = reference_.data(index);
      mask_const_reference mask1 = reference_.mask(index);

      // Create the coordinate system
      vec3<double> m2 = spec_.goniometer().get_rotation_axis();
      vec3<double> s0 = spec_.beam()->get_s0();
      CoordinateSystem cs(m2, s0, s1, phi);

      // Compute the transform
      TransformReverse transform(
          spec_,
          cs,
          sbox.bbox,
          sbox.panel,
          d);

      // Get the transformed shoebox
      data_const_reference p = transform.profile().const_ref();

      // Create the data array
      af::versa< double, af::c_grid<3> > c(sbox.data.accessor());
      std::copy(
          sbox.data.begin(),
          sbox.data.end(),
          c.begin());

      // Create the background array
      af::versa< double, af::c_grid<3> > b(sbox.background.accessor());
      std::copy(
          sbox.background.begin(),
          sbox.background.end(),
          b.begin());

      // Create the mask array
      af::versa< bool, af::c_grid<3> > m(sbox.mask.accessor());

      std::transform(
        sbox.mask.begin(),
        sbox.mask.end(),
        m.begin(),
        detail::check_mask_code(Valid | Foreground));


      // Do the profile fitting
      ProfileFitter<double> fit(
          c.const_ref(),
          b.const_ref(),
          m.const_ref(),
          p,
          1e-3, 100);
      DIALS_ASSERT(fit.niter() < 100);

      // Set the integrated flag
      reflection["intensity.prf.value"] = fit.intensity()[0];
      reflection["intensity.prf.variance"] = fit.variance()[0];
      reflection["intensity.prf.correlation"] = fit.correlation();
      reflection["flags"] = reflection.get<std::size_t>("flags") | af::IntegratedPrf;
    }

    void integrate_detector_space_with_deconvolution(
          af::Reflection &reflection,
          const std::vector<af::Reflection> &adjacent_reflections) const {
      typedef af::const_ref< double, af::c_grid<3> > data_const_reference;
      typedef af::const_ref< bool, af::c_grid<3> > mask_const_reference;





      reflection["flags"] = reflection.get<std::size_t>("flags") & ~af::IntegratedPrf;

      // Get some data
      vec3<double> s1 = reflection.get< vec3<double> >("s1");

      double phi = reflection.get< vec3<double> >("xyzcal.mm")[2];
      vec3<double> xyz = reflection.get< vec3<double> >("xyzcal.px");
      Shoebox<> sbox = reflection.get< Shoebox<> >("shoebox");
      std::size_t flags = reflection.get<std::size_t>("flags");

      bool is_overlapping = false;
      for (std::size_t i =0 ; i < sbox.mask.size(); ++i) {
        if ((sbox.mask[i] & (Valid | Foreground | Overlapped)) == (Valid | Foreground | Overlapped)) {
          is_overlapping = true;
          break;
        }
      }

      /* if (!is_overlapping || adjacent_reflections.size() == 0) { */
      /*   //integrate_detector_space(reflection, adjacent_reflections); */
      /*   return; */
      /* } */

      try {

        DIALS_ASSERT(sbox.is_consistent());


        // Check if we want to use this reflection
        if (check2(flags, sbox) == false) {
          return;
        }

        // Create the data array
        af::versa< double, af::c_grid<3> > c(sbox.data.accessor());
        std::copy(
            sbox.data.begin(),
            sbox.data.end(),
            c.begin());

        // Create the background array
        af::versa< double, af::c_grid<3> > b(sbox.background.accessor());
        std::copy(
            sbox.background.begin(),
            sbox.background.end(),
            b.begin());

        // Create the mask array
        af::versa< bool, af::c_grid<3> > m(sbox.mask.accessor());

        std::transform(
          sbox.mask.begin(),
          sbox.mask.end(),
          m.begin(),
          detail::check_mask_code2(Valid | Foreground));

        // Get the reference profiles
        std::size_t index = sampler_.nearest(sbox.panel, xyz);
        data_const_reference d = reference_.data(index);
        mask_const_reference mask1 = reference_.mask(index);

        // The profile grid
        af::versa< double, af::c_grid<4> > profile(af::c_grid<4>(
              af::tiny<std::size_t,4>(
              adjacent_reflections.size()+1,
              c.accessor()[0],
              c.accessor()[1],
              c.accessor()[2])));

        // Create the coordinate system
        vec3<double> m2 = spec_.goniometer().get_rotation_axis();
        vec3<double> s0 = spec_.beam()->get_s0();
        CoordinateSystem cs(m2, s0, s1, phi);

        // Compute the transform
        TransformReverse transform(
            spec_,
            cs,
            sbox.bbox,
            sbox.panel,
            d);

        // Get the transformed shoebox
        data_const_reference p = transform.profile().const_ref();
        std::copy(p.begin(), p.end(), profile.begin());

        for (std::size_t j = 0; j < adjacent_reflections.size(); ++j) {
          vec3<double> s12 = adjacent_reflections[j].get< vec3<double> >("s1");
          double phi2 = adjacent_reflections[j].get< vec3<double> >("xyzcal.mm")[2];
          CoordinateSystem cs2(m2, s0, s12, phi2);

          // Compute the transform
          TransformReverse transform(
              spec_,
              cs2,
              sbox.bbox,
              sbox.panel,
              d);

          // Get the transformed shoebox
          data_const_reference p2 = transform.profile().const_ref();

          std::copy(p2.begin(), p2.end(), profile.begin() + (j+1)*p2.size());
        }

        /* std::cout << "Num reflections " << adjacent_reflections.size() << std::endl; */
        /* std::cout << "Mask" << std::endl; */
        /* for (std::size_t k = 0; k < m.accessor()[0]; ++k) { */
        /*   for (std::size_t j = 0; j < m.accessor()[1]; ++j) { */
        /*     for (std::size_t i = 0; i < m.accessor()[2]; ++i) { */
        /*       std::cout << m(k,j,i) << ", "; */
        /*     } */
        /*     std::cout << std::endl; */
        /*   } */
        /*   std::cout << std::endl; */
        /*   std::cout << std::endl; */
        /* } */
        /* std::cout << "Data" << std::endl; */
        /* for (std::size_t k = 0; k < m.accessor()[0]; ++k) { */
        /*   for (std::size_t j = 0; j < m.accessor()[1]; ++j) { */
        /*     for (std::size_t i = 0; i < m.accessor()[2]; ++i) { */
        /*       std::cout << c(k,j,i) << ", "; */
        /*     } */
        /*     std::cout << std::endl; */
        /*   } */
        /*   std::cout << std::endl; */
        /*   std::cout << std::endl; */
        /* } */
        /* std::cout << "Background" << std::endl; */
        /* for (std::size_t k = 0; k < m.accessor()[0]; ++k) { */
        /*   for (std::size_t j = 0; j < m.accessor()[1]; ++j) { */
        /*     for (std::size_t i = 0; i < m.accessor()[2]; ++i) { */
        /*       std::cout << b(k,j,i) << ", "; */
        /*     } */
        /*     std::cout << std::endl; */
        /*   } */
        /*   std::cout << std::endl; */
        /*   std::cout << std::endl; */
        /* } */
        /* for (std::size_t l = 0; l < profile.accessor()[0]; ++l) { */
        /*   std::cout << "Profile " << l << std::endl; */
        /*   for (std::size_t k = 0; k < m.accessor()[0]; ++k) { */
        /*     for (std::size_t j = 0; j < m.accessor()[1]; ++j) { */
        /*       for (std::size_t i = 0; i < m.accessor()[2]; ++i) { */
        /*         std::cout << profile(af::tiny<std::size_t,4>(l,k,j,i)) << ", "; */
        /*       } */
        /*       std::cout << std::endl; */
        /*     } */
        /*     std::cout << std::endl; */
        /*     std::cout << std::endl; */
        /*   } */
        /* } */

        // Do the profile fitting
        ProfileFitter<double> fit(
            c.const_ref(),
            b.const_ref(),
            m.const_ref(),
            profile.const_ref(),
            1e-3, 100);
        DIALS_ASSERT(fit.niter() < 100);

        // Set the integrated flag
        reflection["intensity.prf.value"] = fit.intensity()[0];
        reflection["intensity.prf.variance"] = fit.variance()[0];
        reflection["intensity.prf.correlation"] = fit.correlation();
        reflection["flags"] = reflection.get<std::size_t>("flags") | af::IntegratedPrf;

      } catch (dials::error) {
        integrate_detector_space(reflection, adjacent_reflections);
      } catch (std::runtime_error) {
        //integrate_detector_space(reflection, adjacent_reflections);
      }
    }

  protected:

    bool check(const Shoebox<> &sbox) const {

      // Check if the bounding box is in the image
      bool bbox_valid =
        sbox.bbox[0] >= 0 &&
        sbox.bbox[2] >= 0 &&
        sbox.bbox[1] <= spec_.detector()[sbox.panel].get_image_size()[0] &&
        sbox.bbox[3] <= spec_.detector()[sbox.panel].get_image_size()[1];

      // Check if all pixels are valid
      bool pixels_valid = true;
      for (std::size_t i = 0; i < sbox.mask.size(); ++i) {
        if ((sbox.mask[i] & Foreground) &&
            !(sbox.mask[i] & Valid) &&
            !(sbox.mask[i] & Overlapped)) {
          pixels_valid = false;
          break;
        }
      }

      // Return whether to use or not
      return bbox_valid && pixels_valid;
    }

    bool check2(std::size_t flags,
               const Shoebox<> &sbox) const {

      // Check if we want to integrate
      bool integrate = !(flags & af::DontIntegrate);

      // Check if the bounding box is in the image
      bool bbox_valid =
        sbox.bbox[0] >= 0 &&
        sbox.bbox[2] >= 0 &&
        sbox.bbox[1] <= spec_.detector()[sbox.panel].get_image_size()[0] &&
        sbox.bbox[3] <= spec_.detector()[sbox.panel].get_image_size()[1];

      /* // Check if all pixels are valid */
      /* bool pixels_valid = true; */
      /* for (std::size_t i = 0; i < sbox.mask.size(); ++i) { */
      /*   if (sbox.mask[i] & Foreground && !(sbox.mask[i] & Valid)) { */
      /*     pixels_valid = false; */
      /*     break; */
      /*   } */
      /* } */

      // Return whether to use or not
      return integrate && bbox_valid;// && pixels_valid;
    }

    Reference reference_;
    CircleSampler sampler_;
    TransformSpec spec_;
    bool detector_space_;
    bool deconvolution_;

  };

}}

#endif // DIALS_ALGORITHMS_INTEGRATION_ALGORITHMS_H
