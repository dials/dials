#ifndef DIALS_ALGORITHMS_INTEGRATION_PARALLEL_INTEGRATOR_H
#define DIALS_ALGORITHMS_INTEGRATION_PARALLEL_INTEGRATOR_H

/* #include <scitbx/constants.h> */
/* #include <dials/array_family/scitbx_shared_and_versa.h> */
#include <dxtbx/model/beam.h>
#include <dxtbx/model/detector.h>
#include <dxtbx/model/scan.h>
#include <dxtbx/model/crystal.h>
#include <dxtbx/imageset.h>
#include <dials/array_family/reflection_table.h>
#include <dials/array_family/reflection.h>
#include <dials/error.h>
#include <dials/util/thread_pool.h>
#include <dials/algorithms/spot_prediction/pixel_to_miller_index.h>
#include <dials/algorithms/background/glm/robust_poisson_mean.h>
#include <dials/algorithms/integration/sum/summation.h>
#include <dials/algorithms/shoebox/find_overlapping.h>

#include <dials/model/data/mask_code.h>
#include <dials/model/data/shoebox.h>
#include <dials/algorithms/profile_model/gaussian_rs/mask_calculator.h>
#include <dials/algorithms/background/glm/robust_poisson_mean.h>
#include <dials/algorithms/profile_model/modeller/circle_sampler.h>
#include <dials/algorithms/profile_model/gaussian_rs/transform/transform.h>
#include <dials/algorithms/integration/sum/summation.h>
#include <dials/algorithms/integration/fit/fitting.h>
#include <dials/algorithms/centroid/centroid.h>
#include <map>
#include <dials/algorithms/background/glm/creator.h>

namespace dials { namespace algorithms {
  
  using dxtbx::model::Beam;
  using dxtbx::model::Detector;
  using dxtbx::model::Goniometer;
  using dxtbx::model::Scan;

  using dxtbx::ImageSweep;
  using dxtbx::format::Image;
  using dxtbx::format::ImageTile;

  using dials::model::Shoebox;
  using dials::algorithms::profile_model::gaussian_rs::MaskCalculator3D;
  using dials::algorithms::profile_model::gaussian_rs::transform::TransformSpec;
  using dials::algorithms::profile_model::gaussian_rs::CoordinateSystem;
  using dials::algorithms::profile_model::gaussian_rs::transform::TransformReverseNoModel;
  using dials::algorithms::profile_model::gaussian_rs::transform::TransformReverse;
  using dials::algorithms::profile_model::gaussian_rs::transform::TransformForward;
  using dials::algorithms::Creator;

  using dials::model::AdjacencyList;


  class MaskCalculatorIface {
  public:

    virtual void operator()(af::Reflection &reflection, bool adjacent=false) const = 0;

  };

  class BackgroundCalculatorIface {
  public:

    virtual void operator()(af::Reflection &reflection) const = 0;

  };

  class IntensityCalculatorIface {
  public:

    virtual void operator()(af::Reflection &reflection,
                            const std::vector<af::Reflection> &adjacent_reflections) const = 0;

  };


  class MaskCalculator : public MaskCalculatorIface {
  public:

    MaskCalculator(const BeamBase &beam,
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

  class BackgroundCalculator : public BackgroundCalculatorIface {
  public:

    BackgroundCalculator(
          Creator::Model model,
          double tuning_constant,
          std::size_t max_iter,
          std::size_t min_pixels)
      : creator_(
          model,
          tuning_constant,
          max_iter,
          min_pixels) {}

    virtual void operator()(af::Reflection &reflection) const {
      creator_.single(reflection.get< Shoebox<> >("shoebox"));
    }

  protected:

    Creator creator_;
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


  /**
   * A class to store the image data buffer
   */
  class Buffer {
  public:

    /**
     * Initialise the the size of the panels
     * @param detector The detector model
     * @param n The number of images
     */
    Buffer(const Detector &detector, std::size_t n, bool dynamic_mask)
        : dynamic_mask_(dynamic_mask) {
      std::size_t zsize = n;
      DIALS_ASSERT(zsize > 0);
      for (std::size_t i = 0; i < detector.size(); ++i) {
        std::size_t xsize = detector[i].get_image_size()[0];
        std::size_t ysize = detector[i].get_image_size()[1];
        DIALS_ASSERT(xsize > 0);
        DIALS_ASSERT(ysize > 0);
        data_.push_back(
            af::versa< double, af::c_grid<3> >(
              af::c_grid<3>(zsize, ysize, xsize)));
        if (dynamic_mask) {
          mask_.push_back(
            af::versa< bool, af::c_grid<3> >(
              af::c_grid<3>(zsize, ysize, xsize)));
        }
      }
    }

    /**
     * Copy an image to the buffer
     * @param image The image data
     * @param index The image index
     */
    void copy(const Image<double> &data, std::size_t index) {
      DIALS_ASSERT(data.n_tiles() == data_.size());
      for (std::size_t i = 0; i < data.n_tiles(); ++i) {
        copy(data.tile(i).data().const_ref(), data_[i].ref(), index);
      }
    }

    /**
     * Copy an image to the buffer
     * @param image The image data
     * @param index The image index
     */
    void copy(const Image<double> &data, const Image<bool> &mask, std::size_t index) {
      DIALS_ASSERT(dynamic_mask_ == true);
      DIALS_ASSERT(data.n_tiles() == mask.n_tiles());
      DIALS_ASSERT(data.n_tiles() == data_.size());
      DIALS_ASSERT(mask.n_tiles() == mask_.size());
      for (std::size_t i = 0; i < data.n_tiles(); ++i) {
        copy(data.tile(i).data().const_ref(), data_[i].ref(), index);
        copy(mask.tile(i).data().const_ref(), mask_[i].ref(), index);
      }
    }

    /**
     * @param The panel number
     * @returns The buffer for the panel
     */
    af::const_ref< double, af::c_grid<3> > data(std::size_t panel) const {
      DIALS_ASSERT(panel < data_.size());
      return data_[panel].const_ref();
    }

    /**
     * @param The panel number
     * @returns The buffer for the panel
     */
    af::const_ref< bool, af::c_grid<3> > mask(std::size_t panel) const {
      DIALS_ASSERT(dynamic_mask_ == true);
      DIALS_ASSERT(panel < mask_.size());
      return mask_[panel].const_ref();
    }

    /**
     * @returns do we have a dynamic mask?
     */
    bool dynamic_mask() const {
      return dynamic_mask_;
    }

  protected:

    /**
     * Copy the data from 1 panel
     * @param image The image data
     * @param buffer The image buffer
     * @param index The image index
     */
    template <typename T>
    void copy(af::const_ref< T, af::c_grid<2> > src,
              af::ref < T, af::c_grid<3> > dst,
              std::size_t index) {
      std::size_t ysize = src.accessor()[0];
      std::size_t xsize = src.accessor()[1];
      DIALS_ASSERT(index < dst.accessor()[0]);
      DIALS_ASSERT(src.accessor()[0] == dst.accessor()[1]);
      DIALS_ASSERT(src.accessor()[1] == dst.accessor()[2]);
      for (std::size_t j = 0; j < ysize * xsize; ++j) {
        dst[index * (xsize*ysize) + j] = src[j];
      }
    }

    std::vector< af::versa< double, af::c_grid<3> > > data_;
    std::vector< af::versa< bool, af::c_grid<3> > > mask_;
    bool dynamic_mask_;

  };



  /**
   * A class to integrate a single reflection
   */
  class ReflectionIntegrator {
  public:

    /**
     * Initialise the integrator
     * @param compute_mask The mask calculation function
     * @param compute_background The background calculation function
     * @param compute_intensity The intensity calculation function
     * @param buffer The image buffer array
     * @param zstart The first image index
     * @param underload The underload value
     * @param overload The overload value
     */
    ReflectionIntegrator(
          const MaskCalculatorIface &compute_mask,
          const BackgroundCalculatorIface &compute_background,
          const IntensityCalculatorIface &compute_intensity,
          const Buffer &buffer,
          int zstart,
          double underload,
          double overload,
          bool debug)
      : compute_mask_(compute_mask),
        compute_background_(compute_background),
        compute_intensity_(compute_intensity),
        buffer_(buffer),
        zstart_(zstart),
        underload_(underload),
        overload_(overload),
        debug_(debug) {}

    /**
     * Integrate a reflection
     * @param reflection The reflection object
     */
    void operator()(af::Reflection &reflection,
                    std::vector<af::Reflection> &adjacent_reflections) const {

      // Get the panel number
      std::size_t panel = reflection.get<std::size_t>("panel");

      // Extract the shoebox data
      if (buffer_.dynamic_mask()) {
        extract_shoebox(
            buffer_.data(panel),
            buffer_.mask(panel),
            reflection,
            zstart_,
            underload_,
            overload_);
      } else {
        extract_shoebox(
            buffer_.data(panel),
            reflection,
            zstart_,
            underload_,
            overload_);
      }

      // Compute the mask
      compute_mask_(reflection);

      // Set all the bounding boxes of adjacent reflections
      // And compute the mask for these reflections too.
      for (std::size_t i = 0; i < adjacent_reflections.size(); ++i) {
        adjacent_reflections[i]["bbox"] = reflection.get<int6>("bbox");
        adjacent_reflections[i]["shoebox"] = reflection.get< Shoebox<> >("shoebox");
        compute_mask_(adjacent_reflections[i], true);
      }

      // Compute the background
      try {
        compute_background_(reflection);
      } catch (dials::error) {
        return;
      }

      // Compute the centroid
      compute_centroid(reflection);

      // Compute the summed intensity
      compute_summed_intensity(reflection);

      // Compute the profile fitted intensity
      try {
        compute_intensity_(reflection, adjacent_reflections);
      } catch (dials::error) {
        return;
      }

      // Erase the shoebox from the reflection
      if (!debug_) {
        reflection.erase("shoebox");
        for (std::size_t i = 0; i < adjacent_reflections.size(); ++i) {
          adjacent_reflections[i].erase("shoebox");
        }
      }
    }


  protected:

    /**
     * Extract the shoebox
     */
    void extract_shoebox(
          const af::const_ref< double, af::c_grid<3> > &buffer,
          af::Reflection &reflection,
          int zstart,
          double underload,
          double overload) const {
      std::size_t panel = reflection.get<std::size_t>("panel");
      int6 bbox = reflection.get<int6>("bbox");
      Shoebox<> shoebox(panel, bbox);
      shoebox.allocate();
      af::ref< float, af::c_grid<3> > data = shoebox.data.ref();
      af::ref< int,   af::c_grid<3> > mask = shoebox.mask.ref();
      int x0 = bbox[0];
      int x1 = bbox[1];
      int y0 = bbox[2];
      int y1 = bbox[3];
      int z0 = bbox[4];
      int z1 = bbox[5];
      DIALS_ASSERT(x1 > x0);
      DIALS_ASSERT(y1 > y0);
      DIALS_ASSERT(z1 > z0);
      std::size_t zsize = z1 - z0;
      std::size_t ysize = y1 - y0;
      std::size_t xsize = x1 - x0;
      DIALS_ASSERT(zsize == data.accessor()[0]);
      DIALS_ASSERT(ysize == data.accessor()[1]);
      DIALS_ASSERT(xsize == data.accessor()[2]);
      DIALS_ASSERT(shoebox.is_consistent());
      for (std::size_t k = 0; k < zsize; ++k) {
        for (std::size_t j = 0; j < ysize; ++j) {
          for (std::size_t i = 0; i < xsize; ++i) {
            int kk = z0 + k - zstart;
            int jj = y0 + j;
            int ii = x0 + i;
            if (kk >= 0 &&
                jj >= 0 &&
                ii >= 0 &&
                kk < buffer.accessor()[0] &&
                jj < buffer.accessor()[1] &&
                ii < buffer.accessor()[2]) {
              double d = buffer(kk, jj, ii);
              int m = (d > underload && d < overload) ? Valid : 0;
              data(k,j,i) = d;
              mask(k,j,i) = m;
            } else {
              data(k,j,i) = 0;
              mask(k,j,i) = 0;
            }
          }
        }
      }
      reflection["shoebox"] = shoebox;
    }

    /**
     * Extract the shoebox
     */
    void extract_shoebox(
          const af::const_ref< double, af::c_grid<3> > &data_buffer,
          const af::const_ref< bool, af::c_grid<3> > &mask_buffer,
          af::Reflection &reflection,
          int zstart,
          double underload,
          double overload) const {
      DIALS_ASSERT(data_buffer.accessor().all_eq(mask_buffer.accessor()));
      std::size_t panel = reflection.get<std::size_t>("panel");
      int6 bbox = reflection.get<int6>("bbox");
      Shoebox<> shoebox(panel, bbox);
      shoebox.allocate();
      af::ref< float, af::c_grid<3> > data = shoebox.data.ref();
      af::ref< int,   af::c_grid<3> > mask = shoebox.mask.ref();
      int x0 = bbox[0];
      int x1 = bbox[1];
      int y0 = bbox[2];
      int y1 = bbox[3];
      int z0 = bbox[4];
      int z1 = bbox[5];
      DIALS_ASSERT(x1 > x0);
      DIALS_ASSERT(y1 > y0);
      DIALS_ASSERT(z1 > z0);
      std::size_t zsize = z1 - z0;
      std::size_t ysize = y1 - y0;
      std::size_t xsize = x1 - x0;
      DIALS_ASSERT(zsize == data.accessor()[0]);
      DIALS_ASSERT(ysize == data.accessor()[1]);
      DIALS_ASSERT(xsize == data.accessor()[2]);
      DIALS_ASSERT(shoebox.is_consistent());
      for (std::size_t k = 0; k < zsize; ++k) {
        for (std::size_t j = 0; j < ysize; ++j) {
          for (std::size_t i = 0; i < xsize; ++i) {
            int kk = z0 + k - zstart;
            int jj = y0 + j;
            int ii = x0 + i;
            if (kk >= 0 &&
                jj >= 0 &&
                ii >= 0 &&
                kk < data_buffer.accessor()[0] &&
                jj < data_buffer.accessor()[1] &&
                ii < data_buffer.accessor()[2]) {
              double d = data_buffer(kk, jj, ii);
              int m = mask_buffer(kk, jj, ii) ? Valid : 0;
              data(k,j,i) = d;
              mask(k,j,i) = m;
            } else {
              data(k,j,i) = 0;
              mask(k,j,i) = 0;
            }
          }
        }
      }
      reflection["shoebox"] = shoebox;
    }

    /**
     * Compute the centroid
     */
    void compute_centroid(af::Reflection &reflection) const {

      using dials::model::Centroid;

      // Get the shoebox and compute centroid
      Shoebox<> shoebox = reflection.get< Shoebox<> >("shoebox");
      Centroid centroid = shoebox.centroid_foreground_minus_background();

      // Set the centroid values
      reflection["xyzobs.px.value"] = centroid.px.position;
      reflection["xyzobs.px.variance"] = centroid.px.variance;
    }

    /**
     * Compute the summed intensity
     */
    void compute_summed_intensity(af::Reflection &reflection) const {
  
      using dials::model::Intensity;

      // Get flags and reset
      std::size_t flags = reflection.get<std::size_t>("flags");
      flags &= ~af::IntegratedSum;
      flags &= ~af::FailedDuringSummation;

      // Get the shoebox and compute the summed intensity
      Shoebox<> shoebox = reflection.get< Shoebox<> >("shoebox");
      Intensity intensity = shoebox.summed_intensity();

      // Set the intensities
      reflection["intensity.sum.value"] = intensity.observed.value;
      reflection["intensity.sum.variance"] = intensity.observed.variance;
      reflection["background.sum.value"] = intensity.background.value;
      reflection["background.sum.variance"] = intensity.background.variance;

      // Set the appropriate flag
      if (intensity.observed.success) {
        flags |= af::IntegratedSum;
      } else {
        flags |= af::FailedDuringSummation;
      }
      reflection["flags"] = flags;
    }

    const MaskCalculatorIface &compute_mask_;
    const BackgroundCalculatorIface &compute_background_;
    const IntensityCalculatorIface &compute_intensity_;
    const Buffer &buffer_;
    int zstart_;
    double underload_;
    double overload_;
    bool debug_;
  };


  /**
   * a class to sort the indices of all reflections that are fully recorded
   * after a particular image.
   */
  class Lookup {
  public:

    /**
     * @param bbox the bounding boxs
     * @param zstart the first frame number
     * @param n the number of frames
     */
    Lookup(af::const_ref<int6> bbox, int zstart, std::size_t n)
      : indices_(bbox.size()) {

      // fill the index array
      for (std::size_t i = 0; i < indices_.size(); ++i) {
        indices_[i] = i;
      }

      // sort the indices by the final frame number
      std::sort(indices_.begin(), indices_.end(), sort_by_frame(bbox));
      DIALS_ASSERT(bbox[indices_.front()][5] - zstart >= 1);
      DIALS_ASSERT(bbox[indices_.back()][5]  - zstart <= n);

      // create an offset array that records the positions in the index array
      // where the frame increments such that
      // offset[i], offset[i+1] gives the range of indices for a frame.
      std::size_t i = 0;
      offset_.push_back(0);
      for (std::size_t j = 0; j < n; ++j) {
        while (i < indices_.size() && bbox[indices_[i]][5] - zstart <= j+1) i++;
        offset_.push_back(i);
      }
      DIALS_ASSERT(offset_.size() == n + 1);
      DIALS_ASSERT(offset_.back() == indices_.size());
    }

    /**
     * @param z the frame number
     * @returns the indices for a given frame
     */
    af::const_ref<std::size_t> indices(std::size_t z) const {
      DIALS_ASSERT(z < offset_.size()-1);
      DIALS_ASSERT(offset_[z+1] >= offset_[z]);
      std::size_t i = offset_[z];
      std::size_t n = offset_[z+1] - i;
      return af::const_ref<std::size_t>(&indices_[i], n);
    }

  private:

    /**
     * helper function to sort by final bbox frame
     */
    struct sort_by_frame {
      af::const_ref<int6> bbox_;
      sort_by_frame(af::const_ref<int6> bbox) : bbox_(bbox) {}
      bool operator()(std::size_t a, std::size_t b) const {
        return bbox_[a][5] < bbox_[b][5];
      }
    };

    std::vector<std::size_t> indices_;
    std::vector<std::size_t> offset_;
  };



  /**
   * A class to perform parallel integration
   */
  class Integrator {
  public:

    /**
     * Do the integration
     * @param reflections The reflection table
     * @param imageset The imageset
     * @param compute_mask The mask calulcation function
     * @param compute_background The background calculation function
     * @param compute_intensity The intensity calculation function
     * @param nthreads The number of parallel threads
     */
    Integrator(
          af::reflection_table reflections,
          ImageSweep imageset,
          const MaskCalculatorIface &compute_mask,
          const BackgroundCalculatorIface &compute_background,
          const IntensityCalculatorIface &compute_intensity,
          std::size_t nthreads,
          bool debug) {
  
      using dials::algorithms::shoebox::find_overlapping_multi_panel;

      // Check the input
      DIALS_ASSERT(nthreads > 0);

      // Check the models
      DIALS_ASSERT(imageset.get_detector() != NULL);
      DIALS_ASSERT(imageset.get_scan() != NULL);
      Detector detector = *imageset.get_detector();
      Scan scan = *imageset.get_scan();

      // Get the size of the data buffer needed
      std::size_t zsize = imageset.size();
      DIALS_ASSERT(zsize > 0);

      // Get the starting frame and the underload/overload values
      int zstart = scan.get_array_range()[0];
      double underload = detector[0].get_trusted_range()[0];
      double overload  = detector[0].get_trusted_range()[1];
      DIALS_ASSERT(underload < overload);
      for (std::size_t i = 1; i < detector.size(); ++i) {
        DIALS_ASSERT(underload == detector[i].get_trusted_range()[0]);
        DIALS_ASSERT(overload  == detector[i].get_trusted_range()[1]);
      }

      // Get the reflection flags and bbox
      af::const_ref<std::size_t> panel = reflections.get<std::size_t>("panel").const_ref();
      af::const_ref<int6> bbox = reflections.get<int6>("bbox").const_ref();
      af::ref<std::size_t> flags = reflections.get<std::size_t>("flags").ref();

      // Find the overlapping reflections
      AdjacencyList overlaps = find_overlapping_multi_panel(bbox, panel);

      // Allocate the array for the image data
      Buffer buffer(detector, zsize, false);

      // Transform reflection data from column major to row major. The
      // reason for doing this is because we want to process each reflection
      // in parallel. Therefore passing a reflection object is better than
      // passing the whole table. Additionally, because the reflection table
      // uses a std::map, it may not be thread safe, so we want to avoid
      // accessing this across multiple threads.
      af::shared<af::Reflection> reflection_array =
        reflection_table_to_array(reflections);

      // The lookup class gives the indices of reflections whose bounding boxes
      // are complete at a given image. This is used to submit reflections for
      // integration after each image is processed.
      Lookup lookup(bbox, zstart, zsize);

      // Create the reflection integrator. This class is called for each
      // reflection to integrate the data
      ReflectionIntegrator integrator(
          compute_mask,
          compute_background,
          compute_intensity,
          buffer,
          zstart,
          underload,
          overload,
          debug);

      // Do the integration
      process(
          lookup,
          integrator,
          buffer,
          reflection_array.ref(),
          overlaps,
          imageset,
          bbox,
          flags,
          nthreads);

      // Transform the row major reflection array to the reflection table
      reflections_ = reflection_table_from_array(reflection_array.const_ref());
    }

    /**
     * @returns The integrated reflections
     */
    af::reflection_table reflections() const {
      return reflections_;
    }

  protected:

    /**
     * Reset the reflection flags
     */
    void reset_flags(af::ref<std::size_t> flags) const {
      for (std::size_t i = 0; i < flags.size(); ++i) {
        flags[i] &= ~af::IntegratedSum;
        flags[i] &= ~af::IntegratedPrf;
      }
    }

    /**
     * Do the processing
     */
    void process(
        const Lookup &lookup,
        const ReflectionIntegrator &integrator,
        Buffer &buffer,
        af::ref< af::Reflection > reflections,
        const AdjacencyList &overlaps,
        ImageSweep imageset,
        af::const_ref<int6> bbox,
        af::const_ref<std::size_t> flags,
        std::size_t nthreads) const {
  
      using dials::util::ThreadPool;

      // Construct a list of adjacent reflections for each given reflection
      std::vector< std::vector<af::Reflection> > adjacent_reflections(reflections.size());
      for (std::size_t i = 0; i < reflections.size(); ++i) {
        if ((flags[i] & af::DontIntegrate) == 0) {
          adjacent_reflections[i].reserve(overlaps.vertex_num_edges(i));
          AdjacencyList::edge_iterator_range edges = overlaps.edges(i);
          for (AdjacencyList::edge_iterator it = edges.first; it != edges.second; ++it) {
            DIALS_ASSERT(it->first == i);
            DIALS_ASSERT(it->second < reflections.size());
            adjacent_reflections[i].push_back(reflections[it->second]);
          }
        }
      }

      // Create the thread pool
      ThreadPool pool(nthreads);

      // Get the size of the array
      int zstart = imageset.get_scan()->get_array_range()[0];
      std::size_t zsize = imageset.size();

      // Loop through all the images
      for (std::size_t i = 0; i < zsize; ++i) {
        std::cout << zstart + i << std::endl;

        // Copy the image to the buffer
        if (buffer.dynamic_mask()) {
          buffer.copy(imageset.get_corrected_data(i), imageset.get_mask(i), i);
        } else {
          buffer.copy(imageset.get_corrected_data(i), i);
        }

        // Get the reflections recorded at this point
        af::const_ref<std::size_t> indices = lookup.indices(i);

        // Iterate through the reflection indices
        for (std::size_t j = 0; j < indices.size(); ++j) {

          // Get the reflection index
          std::size_t k = indices[j];

          // Check that the reflection bounding box is within the range of
          // images that have been read
          DIALS_ASSERT(bbox[k][5] <= zstart + i + 1);

          // Ignore if we're not integrating this reflection
          if (flags[k] & af::DontIntegrate) {
            continue;
          }

          // Post the integration job
          pool.post(
            boost::bind(
              &ReflectionIntegrator::operator(),
              integrator,
              boost::ref(reflections[k]),
              boost::ref(adjacent_reflections[k])));
        }
      }

      // Wait for all the integration jobs to complete
      pool.wait();
    }

    af::reflection_table reflections_;
  };

}}

#endif // DIALS_ALGORITHMS_INTEGRATION_PARALLEL_INTEGRATOR_H
