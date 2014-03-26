/*
 * transform.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_REFLEXION_BASIS_MAP_PIXELS2_H
#define DIALS_ALGORITHMS_REFLEXION_BASIS_MAP_PIXELS2_H

#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include <scitbx/array_family/tiny_types.h>
#include <dxtbx/model/beam.h>
#include <dxtbx/model/detector.h>
#include <dxtbx/model/goniometer.h>
#include <dxtbx/model/scan.h>
#include <dials/algorithms/polygon/spatial_interpolation.h>
#include <dials/algorithms/reflection_basis/coordinate_system.h>
#include <dials/algorithms/reflection_basis/map_frames.h>
#include <dials/algorithms/reflection_basis/beam_vector_map.h>
#include <dials/model/data/shoebox.h>

namespace dials { namespace algorithms { namespace reflection_basis {
    namespace transform {

  using scitbx::vec2;
  using scitbx::vec3;
  using scitbx::af::int2;
  using scitbx::af::int3;
  using scitbx::af::double3;
  using scitbx::af::int6;
  using dxtbx::model::Beam;
  using dxtbx::model::Detector;
  using dxtbx::model::Goniometer;
  using dxtbx::model::Scan;
  using dials::model::Foreground;
  using dials::model::Valid;
  using dials::model::Shoebox;
  using dials::algorithms::polygon::spatial_interpolation::vert4;
  using dials::algorithms::polygon::spatial_interpolation::Match;
  using dials::algorithms::polygon::spatial_interpolation::quad_to_grid;

  /**
   * A class to construct the specification for the transform. Once instantiated
   * this object can be reused to transform lots of reflections.
   */
  template <typename FloatType = double>
  class TransformSpec {
  public:

    typedef FloatType float_type;
    typedef MapFramesForward<FloatType> map_frames_type;

    /**
     * Initialise the class
     * @param beam The beam model
     * @param detector The detector model
     * @param gonio The goniometer model
     * @param scan The scan model
     * @param sigma_b The beam divergence
     * @param sigma_m The crystal mosaicity
     * @param n_sigma The number of standard deviations
     * @param grid_size The size of the reflection basis grid
     */
    TransformSpec(const Beam &beam, const Detector &detector,
                  const Goniometer &gonio, const Scan &scan,
                  double sigma_b, double sigma_m, double n_sigma,
                  std::size_t grid_size)
      : n_div_(5),
        s0_(beam.get_s0()),
        m2_(gonio.get_rotation_axis().normalize()),
        image_size_(detector[0].get_image_size()[1],
                    detector[0].get_image_size()[0]),
        grid_size_(2*grid_size+1, 2*grid_size+1, 2*grid_size+1),
        step_size_(sigma_m * n_sigma / (grid_size + 0.5),
                   sigma_b * n_sigma / (grid_size + 0.5),
                   sigma_b * n_sigma / (grid_size + 0.5)),
        grid_centre_(grid_size + 0.5, grid_size + 0.5, grid_size + 0.5),
        s1_map_(beam_vector_map(detector, beam, n_div_, false)),
        map_frames_(scan.get_oscillation()[0],
                    scan.get_oscillation()[1],
                    sigma_m, n_sigma, grid_size) {
      DIALS_ASSERT(detector.size() == 1);
      DIALS_ASSERT(image_size_.all_gt(0));
      DIALS_ASSERT(step_size_.all_gt(0));
      DIALS_ASSERT(grid_size_.all_gt(0));
    }

    /** @returns The number of pixel sub division */
    std::size_t n_div() const {
      return n_div_;
    }

    /** @ returns the rotation angle */
    vec3<double> m2() const {
      return m2_;
    }

    /** @returns the incident beam vector */
    vec3<double> s0() const {
      return s0_;
    }

    /** @returns the image size */
    int2 image_size() const {
      return image_size_;
    }

    /** @returns the grid size */
    int3 grid_size() const {
      return grid_size_;
    }

    /** @returns the grid step size */
    double3 step_size() const {
      return step_size_;
    }

    /** @returns the grid centre */
    double3 grid_centre() const {
      return grid_centre_;
    }

    /** @returns the beam vector lookup map */
    af::versa< vec3<double>, af::c_grid<2> > s1_map() const {
      return s1_map_;
    }

    /** @returns the frame mapping fraction array */
    af::versa< FloatType, af::c_grid<2> > map_frames(
        vec2<int> frames, double phi, double zeta) const {
      return map_frames_(frames, phi, zeta);
    }

  private:
    std::size_t n_div_;
    vec3<double> s0_;
    vec3<double> m2_;
    int2 image_size_;
    int3 grid_size_;
    double3 step_size_;
    double3 grid_centre_;
    af::versa< vec3<double>, af::c_grid<2> > s1_map_;
    MapFramesForward<FloatType> map_frames_;
  };


  /**
   * A class to perform the local coordinate transform for a single reflection.
   * The class has a number of different constructors to allow the transform
   * to be done with lots of different inputs.
   *
   * Example:
   *
   *  from dials.algorithms.reflection_basis import transform
   *  forward = transform.Forward(spec, reflection)
   *  print forward.profile()
   *  print forward.background()
   */
  template <typename FloatType = double>
  class Forward {
  public:

    typedef FloatType float_type;
    typedef TransformSpec<FloatType> transform_spec_type;

    Forward() {}

    Forward(const TransformSpec<FloatType> &spec,
            const vec3<double> &s1, double phi, int6 bbox,
            const af::const_ref< FloatType, af::c_grid<3> > &image,
            const af::const_ref< bool, af::c_grid<3> > &mask) {
      init(spec, s1, phi, bbox);
      call(image, mask);
    }

    Forward(const TransformSpec<FloatType> &spec,
            const vec3<double> &s1, double phi, int6 bbox,
            const af::const_ref< FloatType, af::c_grid<3> > &image,
            const af::const_ref< FloatType, af::c_grid<3> > &bkgrd,
            const af::const_ref< bool, af::c_grid<3> > &mask) {
      init(spec, s1, phi, bbox);
      call(image, bkgrd, mask);
    }

    Forward(const TransformSpec<FloatType> &spec,
            const vec3<double> &s1, double phi,
            const Shoebox<FloatType> &shoebox) {
      init(spec, s1, phi, shoebox.bbox);
      call(shoebox);
    }

    Forward(const TransformSpec<FloatType> &spec,
            const CoordinateSystem &cs, int6 bbox,
            const af::const_ref< FloatType, af::c_grid<3> > &image,
            const af::const_ref< bool, af::c_grid<3> > &mask) {
      init(spec, cs, bbox);
      call(image, mask);
    }

    Forward(const TransformSpec<FloatType> &spec,
            const CoordinateSystem &cs, int6 bbox,
            const af::const_ref< FloatType, af::c_grid<3> > &image,
            const af::const_ref< FloatType, af::c_grid<3> > &bkgrd,
            const af::const_ref< bool, af::c_grid<3> > &mask) {
      init(spec, cs, bbox);
      call(image, bkgrd, mask);
    }

    Forward(const TransformSpec<FloatType> &spec,
            const CoordinateSystem &cs,
            const Shoebox<FloatType> &shoebox) {
      init(spec, cs, shoebox.bbox);
      call(shoebox);
    }

    /** @returns The transformed profile */
    af::versa< FloatType, af::c_grid<3> > profile() const {
      return profile_;
    }

    /** @returns The transformed background (if set) */
    af::versa< FloatType, af::c_grid<3> > background() const {
      return background_;
    }

    /** @returns The z fraction */
    af::versa< FloatType, af::c_grid<2> > zfraction() const {
      return zfraction_arr_;
    }

  private:

    /** Initialise using the beam vector and rotation angle */
    void init(const TransformSpec<FloatType> &spec,
              const vec3<double> &s1, double phi, int6 bbox) {
      CoordinateSystem cs(spec.m2(), spec.s0(), s1, phi);
      init(spec, cs, bbox);
    }

    /** Initialise using a coordinate system struct */
    void init(const TransformSpec<FloatType> &spec,
              const CoordinateSystem &cs, int6 bbox) {

      // Initialise some stuff
      x0_ = bbox[0];
      y0_ = bbox[2];
      shoebox_size_ = int3(bbox[5]-bbox[4], bbox[3]-bbox[2], bbox[1]-bbox[0]);
      DIALS_ASSERT(shoebox_size_.all_gt(0));
      DIALS_ASSERT(bbox[0] >= 0 && bbox[2] >= 0);
      DIALS_ASSERT(bbox[1] <= spec.image_size()[1]);
      DIALS_ASSERT(bbox[3] <= spec.image_size()[0]);
      step_size_ = spec.step_size();
      grid_size_ = spec.grid_size();
      grid_cent_ = spec.grid_centre();
      s1_ = cs.s1();
      DIALS_ASSERT(s1_.length() > 0);
      e1_ = cs.e1_axis() / s1_.length();
      e2_ = cs.e2_axis() / s1_.length();
      s1_map_arr_ = spec.s1_map();
      s1_map_ = s1_map_arr_.const_ref();
      n_div_ = spec.n_div();

      // Calculate the fraction of intensity contributed from each data
      // frame to each grid coordinate
      vec2<int> zrange(bbox[4], bbox[5]);
      zfraction_arr_ = spec.map_frames(zrange, cs.phi(), cs.zeta());
      zfraction_ = zfraction_arr_.const_ref();
    }

    /**
     * Map the pixel values from the input image to the output grid.
     * @param image The image to transform
     * @param mask The mask accompanying the image
     */
    void call(const af::const_ref< FloatType, af::c_grid<3> > &image,
              const af::const_ref< bool, af::c_grid<3> > &mask) {

      // Check the input
      DIALS_ASSERT(image.accessor().all_eq(shoebox_size_));
      DIALS_ASSERT(image.accessor().all_eq(mask.accessor()));

      // Initialise the profile arrays
      af::c_grid<3> accessor(grid_size_);
      profile_ = af::versa< FloatType, af::c_grid<3> >(accessor, 0.0);

      // Loop through all the points in the shoebox. Calculate the polygon
      // formed by the pixel in the local coordinate system. Find the points
      // on the grid which intersect with the polygon and the fraction of the
      // pixel area shared with each grid point. For each intersection, loop
      // through the frames, mapping the fraction of the pixel value in each
      // frame to the grid point.
      std::size_t ndiv = n_div_;
      double fraction = 1.0 / (ndiv * ndiv);
      for (std::size_t j = 0; j < shoebox_size_[1]; ++j) {
        for (std::size_t i = 0; i < shoebox_size_[2]; ++i) {
          for (std::size_t jj = 0; jj < ndiv; ++jj) {
            for (std::size_t ii = 0; ii < ndiv; ++ii) {
              std::size_t yy = (y0_ + j) * ndiv + jj;
              std::size_t xx = (x0_ + i) * ndiv + ii;
              vec2<double> gxy = gc(yy, xx);
              int gj = (int)gxy[0];
              int gi = (int)gxy[1];
              if (gj < 0 || gj >= grid_size_[1] ||
                  gi < 0 || gi >= grid_size_[2]) {
                continue;
              }
              for (std::size_t k = 0; k < shoebox_size_[0]; ++k) {
                if (mask(k,j,i)) {
                  FloatType value = image(k,j,i) * fraction;
                  for (std::size_t gk = 0; gk < grid_size_[0]; ++gk) {
                    profile_(gk, gj, gi) += value * zfraction_(k, gk);
                  }
                }
              }
            }
          }
        }
      }
    }

    /**
     * Map the pixel values from the input image to the output grid.
     * @param image The image to transform
     * @param bkgrd The background image to transform
     * @param mask The mask accompanying the image
     */
    void call(const af::const_ref< FloatType, af::c_grid<3> > &image,
              const af::const_ref< FloatType, af::c_grid<3> > &bkgrd,
              const af::const_ref< bool, af::c_grid<3> > &mask) {

      // Check the input
      DIALS_ASSERT(image.accessor().all_eq(shoebox_size_));
      DIALS_ASSERT(image.accessor().all_eq(mask.accessor()));
      DIALS_ASSERT(image.accessor().all_eq(bkgrd.accessor()));

      // Initialise the profile arrays
      af::c_grid<3> accessor(grid_size_);
      profile_ = af::versa< FloatType, af::c_grid<3> >(accessor, 0.0);
      background_ = af::versa< FloatType, af::c_grid<3> >(accessor, 0.0);

      // Loop through all the points in the shoebox. Calculate the polygon
      // formed by the pixel in the local coordinate system. Find the points
      // on the grid which intersect with the polygon and the fraction of the
      // pixel area shared with each grid point. For each intersection, loop
      // through the frames, mapping the fraction of the pixel value in each
      // frame to the grid point.
      std::size_t ndiv = n_div_;
      double fraction = 1.0 / (ndiv * ndiv);
      for (std::size_t j = 0; j < shoebox_size_[1]; ++j) {
        for (std::size_t i = 0; i < shoebox_size_[2]; ++i) {
          for (std::size_t jj = 0; jj < ndiv; ++jj) {
            for (std::size_t ii = 0; ii < ndiv; ++ii) {
              std::size_t yy = (y0_ + j) * ndiv + jj;
              std::size_t xx = (x0_ + i) * ndiv + ii;
              vec2<double> gxy = gc(yy, xx);
              int gj = (int)gxy[0];
              int gi = (int)gxy[1];
              if (gj < 0 || gj >= grid_size_[1] ||
                  gi < 0 || gi >= grid_size_[2]) {
                continue;
              }
              for (std::size_t k = 0; k < shoebox_size_[0]; ++k) {
                if (mask(k,j,i)) {
                  FloatType ivalue = image(k, j, i) * fraction;
                  FloatType bvalue = bkgrd(k, j, i) * fraction;
                  for (std::size_t gk = 0; gk < grid_size_[0]; ++gk) {
                    FloatType zf = zfraction_(k, gk);
                    profile_(gk, gj, gi) += ivalue * zf;
                    background_(gk, gj, gi) += bvalue * zf;
                  }
                }
              }
            }
          }
        }
      }
    }

    /**
     * Call the transform with the shoebox
     */
    void call(const Shoebox<> &shoebox) {
      af::versa< bool, af::c_grid<3> > mask(shoebox.mask.accessor());
      af::ref< bool, af::c_grid<3> > mask_ref = mask.ref();
      af::const_ref< int, af::c_grid<3> > temp_ref = shoebox.mask.const_ref();
      for (std::size_t i = 0; i < mask_ref.size(); ++i) {
        mask_ref[i] = (temp_ref[i] & Valid && temp_ref[i] & Foreground);
      }
      call(shoebox.data.const_ref(), shoebox.background.const_ref(), mask_ref);
    }

    /**
     * Get a grid coordinate from an image coordinate
     * @param j The y index
     * @param i The x index
     * @returns The grid (c1, c2) index
     */
    vec2<double> gc(std::size_t j, std::size_t i) const {
      DIALS_ASSERT(j < s1_map_.accessor()[0] && i < s1_map_.accessor()[1]);
      vec3<double> ds = s1_map_(j, i) - s1_;
      return vec2<double>(grid_cent_[2] + (e1_ * ds) / step_size_[2],
                          grid_cent_[1] + (e2_ * ds) / step_size_[1]);
    }

    std::size_t n_div_;
    int x0_, y0_;
    int3 shoebox_size_;
    int3 grid_size_;
    double3 step_size_;
    double3 grid_cent_;
    vec3<double> s1_, e1_, e2_;
    af::versa< vec3<double>, af::c_grid<2> > s1_map_arr_;
    af::const_ref< vec3<double>, af::c_grid<2> > s1_map_;
    af::versa< FloatType, af::c_grid<3> > profile_;
    af::versa< FloatType, af::c_grid<3> > background_;
    af::versa< FloatType, af::c_grid<2> > zfraction_arr_;
    af::const_ref< FloatType, af::c_grid<2> > zfraction_;
  };

  /**
   * A class to construct the specification for the transform. Once instantiated
   * this object can be reused to transform lots of reflections.
   */
  //template <typename FloatType = double>
  //class TransformSpec {
  //public:

    //typedef FloatType float_type;
    //typedef MapFramesForward<FloatType> map_frames_type;

    /**
     * Initialise the class
     * @param beam The beam model
     * @param detector The detector model
     * @param gonio The goniometer model
     * @param scan The scan model
     * @param sigma_b The beam divergence
     * @param sigma_m The crystal mosaicity
     * @param n_sigma The number of standard deviations
     * @param grid_size The size of the reflection basis grid
     */
    //TransformSpec(const Beam &beam, const Detector &detector,
                  //const Goniometer &gonio, const Scan &scan,
                  //double sigma_b, double sigma_m, double n_sigma,
                  //std::size_t grid_size)
      //: s0_(beam.get_s0()),
        //m2_(gonio.get_rotation_axis().normalize()),
        //image_size_(detector[0].get_image_size()[1],
                    //detector[0].get_image_size()[0]),
        //grid_size_(2*grid_size+1, 2*grid_size+1, 2*grid_size+1),
        //step_size_(sigma_m * n_sigma / (grid_size + 0.5),
                   //sigma_b * n_sigma / (grid_size + 0.5),
                   //sigma_b * n_sigma / (grid_size + 0.5)),
        //grid_centre_(grid_size + 0.5, grid_size + 0.5, grid_size + 0.5),
        //s1_map_(beam_vector_map(detector, beam, true)),
        //map_frames_(scan.get_oscillation()[0],
                    //scan.get_oscillation()[1],
                    //sigma_m, n_sigma, grid_size) {
      //DIALS_ASSERT(detector.size() == 1);
      //DIALS_ASSERT(image_size_.all_gt(0));
      //DIALS_ASSERT(step_size_.all_gt(0));
      //DIALS_ASSERT(grid_size_.all_gt(0));
    //}

    //[>* @ returns the rotation angle <]
    //vec3<double> m2() const {
      //return m2_;
    //}

    //[>* @returns the incident beam vector <]
    //vec3<double> s0() const {
      //return s0_;
    //}

    //[>* @returns the image size <]
    //int2 image_size() const {
      //return image_size_;
    //}

    //[>* @returns the grid size <]
    //int3 grid_size() const {
      //return grid_size_;
    //}

    //[>* @returns the grid step size <]
    //double3 step_size() const {
      //return step_size_;
    //}

    //[>* @returns the grid centre <]
    //double3 grid_centre() const {
      //return grid_centre_;
    //}

    //[>* @returns the beam vector lookup map <]
    //af::versa< vec3<double>, af::c_grid<2> > s1_map() const {
      //return s1_map_;
    //}

    //[>* @returns the frame mapping fraction array <]
    //af::versa< FloatType, af::c_grid<2> > map_frames(
        //vec2<int> frames, double phi, double zeta) const {
      //return map_frames_(frames, phi, zeta);
    //}

  //private:
    //vec3<double> s0_;
    //vec3<double> m2_;
    //int2 image_size_;
    //int3 grid_size_;
    //double3 step_size_;
    //double3 grid_centre_;
    //af::versa< vec3<double>, af::c_grid<2> > s1_map_;
    //MapFramesForward<FloatType> map_frames_;
  //};


  /**
   * A class to perform the local coordinate transform for a single reflection.
   * The class has a number of different constructors to allow the transform
   * to be done with lots of different inputs.
   *
   * Example:
   *
   *  from dials.algorithms.reflection_basis import transform
   *  forward = transform.Forward(spec, reflection)
   *  print forward.profile()
   *  print forward.background()
   */
  //template <typename FloatType = double>
  //class Forward {
  //public:

    //typedef FloatType float_type;
    //typedef TransformSpec<FloatType> transform_spec_type;

    //Forward() {}

    //Forward(const TransformSpec<FloatType> &spec,
            //const vec3<double> &s1, double phi, int6 bbox,
            //const af::const_ref< FloatType, af::c_grid<3> > &image,
            //const af::const_ref< bool, af::c_grid<3> > &mask) {
      //init(spec, s1, phi, bbox);
      //call(image, mask);
    //}

    //Forward(const TransformSpec<FloatType> &spec,
            //const vec3<double> &s1, double phi, int6 bbox,
            //const af::const_ref< FloatType, af::c_grid<3> > &image,
            //const af::const_ref< FloatType, af::c_grid<3> > &bkgrd,
            //const af::const_ref< bool, af::c_grid<3> > &mask) {
      //init(spec, s1, phi, bbox);
      //call(image, bkgrd, mask);
    //}

    //Forward(const TransformSpec<FloatType> &spec,
            //const vec3<double> &s1, double phi,
            //const Shoebox<FloatType> &shoebox) {
      //init(spec, s1, phi, shoebox.bbox);
      //call(shoebox);
    //}

    //Forward(const TransformSpec<FloatType> &spec,
            //const CoordinateSystem &cs, int6 bbox,
            //const af::const_ref< FloatType, af::c_grid<3> > &image,
            //const af::const_ref< bool, af::c_grid<3> > &mask) {
      //init(spec, cs, bbox);
      //call(image, mask);
    //}

    //Forward(const TransformSpec<FloatType> &spec,
            //const CoordinateSystem &cs, int6 bbox,
            //const af::const_ref< FloatType, af::c_grid<3> > &image,
            //const af::const_ref< FloatType, af::c_grid<3> > &bkgrd,
            //const af::const_ref< bool, af::c_grid<3> > &mask) {
      //init(spec, cs, bbox);
      //call(image, bkgrd, mask);
    //}

    //Forward(const TransformSpec<FloatType> &spec,
            //const CoordinateSystem &cs,
            //const Shoebox<FloatType> &shoebox) {
      //init(spec, cs, shoebox.bbox);
      //call(shoebox);
    //}

    //[>* @returns The transformed profile <]
    //af::versa< FloatType, af::c_grid<3> > profile() const {
      //return profile_;
    //}

    //[>* @returns The transformed background (if set) <]
    //af::versa< FloatType, af::c_grid<3> > background() const {
      //return background_;
    //}

    //[>* @returns The z fraction <]
    //af::versa< FloatType, af::c_grid<2> > zfraction() const {
      //return zfraction_arr_;
    //}

  //private:

    //[>* Initialise using the beam vector and rotation angle <]
    //void init(const TransformSpec<FloatType> &spec,
              //const vec3<double> &s1, double phi, int6 bbox) {
      //CoordinateSystem cs(spec.m2(), spec.s0(), s1, phi);
      //init(spec, cs, bbox);
    //}

    //[>* Initialise using a coordinate system struct <]
    //void init(const TransformSpec<FloatType> &spec,
              //const CoordinateSystem &cs, int6 bbox) {

      //// Initialise some stuff
      //x0_ = bbox[0];
      //y0_ = bbox[2];
      //shoebox_size_ = int3(bbox[5]-bbox[4], bbox[3]-bbox[2], bbox[1]-bbox[0]);
      //DIALS_ASSERT(shoebox_size_.all_gt(0));
      //DIALS_ASSERT(bbox[0] >= 0 && bbox[2] >= 0);
      //DIALS_ASSERT(bbox[1] <= spec.image_size()[1]);
      //DIALS_ASSERT(bbox[3] <= spec.image_size()[0]);
      //step_size_ = spec.step_size();
      //grid_size_ = spec.grid_size();
      //grid_cent_ = spec.grid_centre();
      //s1_ = cs.s1();
      //DIALS_ASSERT(s1_.length() > 0);
      //e1_ = cs.e1_axis() / s1_.length();
      //e2_ = cs.e2_axis() / s1_.length();
      //s1_map_arr_ = spec.s1_map();
      //s1_map_ = s1_map_arr_.const_ref();

      //// Calculate the fraction of intensity contributed from each data
      //// frame to each grid coordinate
      //vec2<int> zrange(bbox[4], bbox[5]);
      //zfraction_arr_ = spec.map_frames(zrange, cs.phi(), cs.zeta());
      //zfraction_ = zfraction_arr_.const_ref();
    //}

    /**
     * Map the pixel values from the input image to the output grid.
     * @param image The image to transform
     * @param mask The mask accompanying the image
     */
    //void call(const af::const_ref< FloatType, af::c_grid<3> > &image,
              //const af::const_ref< bool, af::c_grid<3> > &mask) {

      //// Check the input
      //DIALS_ASSERT(image.accessor().all_eq(shoebox_size_));
      //DIALS_ASSERT(image.accessor().all_eq(mask.accessor()));

      //// Initialise the profile arrays
      //af::c_grid<3> accessor(grid_size_);
      //profile_ = af::versa< FloatType, af::c_grid<3> >(accessor, 0.0);

      //// Loop through all the points in the shoebox. Calculate the polygon
      //// formed by the pixel in the local coordinate system. Find the points
      //// on the grid which intersect with the polygon and the fraction of the
      //// pixel area shared with each grid point. For each intersection, loop
      //// through the frames, mapping the fraction of the pixel value in each
      //// frame to the grid point.
      //int2 grid_size2(grid_size_[1], grid_size_[2]);
      //for (std::size_t j = 0; j < shoebox_size_[1]; ++j) {
        //for (std::size_t i = 0; i < shoebox_size_[2]; ++i) {
          //vert4 input(gc(j, i), gc(j, i+1), gc(j+1, i+1), gc(j+1, i));
          //af::shared<Match> matches = quad_to_grid(input, grid_size2, 0);
          //for (int m = 0; m < matches.size(); ++m) {
            //FloatType fraction = matches[m].fraction;
            //int index = matches[m].out;
            //int ii = index % grid_size_[2];
            //int jj = index / grid_size_[2];
            //for (int k = 0; k < shoebox_size_[0]; ++k) {
              //if (mask(k, j, i)) {
                //FloatType value = image(k, j, i) * fraction;
                //for (int kk = 0; kk < grid_size_[0]; ++kk) {
                  //profile_(kk, jj, ii) += value * zfraction_(k, kk);
                //}
              //}
            //}
          //}
        //}
      //}
    //}

    /**
     * Map the pixel values from the input image to the output grid.
     * @param image The image to transform
     * @param bkgrd The background image to transform
     * @param mask The mask accompanying the image
     */
    //void call(const af::const_ref< FloatType, af::c_grid<3> > &image,
              //const af::const_ref< FloatType, af::c_grid<3> > &bkgrd,
              //const af::const_ref< bool, af::c_grid<3> > &mask) {

      //// Check the input
      //DIALS_ASSERT(image.accessor().all_eq(shoebox_size_));
      //DIALS_ASSERT(image.accessor().all_eq(mask.accessor()));
      //DIALS_ASSERT(image.accessor().all_eq(bkgrd.accessor()));

      //// Initialise the profile arrays
      //af::c_grid<3> accessor(grid_size_);
      //profile_ = af::versa< FloatType, af::c_grid<3> >(accessor, 0.0);
      //background_ = af::versa< FloatType, af::c_grid<3> >(accessor, 0.0);

      //// Loop through all the points in the shoebox. Calculate the polygon
      //// formed by the pixel in the local coordinate system. Find the points
      //// on the grid which intersect with the polygon and the fraction of the
      //// pixel area shared with each grid point. For each intersection, loop
      //// through the frames, mapping the fraction of the pixel value in each
      //// frame to the grid point.
      //int2 grid_size2(grid_size_[1], grid_size_[2]);
      //for (std::size_t j = 0; j < shoebox_size_[1]; ++j) {
        //for (std::size_t i = 0; i < shoebox_size_[2]; ++i) {
          //vert4 input(gc(j, i), gc(j, i+1), gc(j+1, i+1), gc(j+1, i));
          //af::shared<Match> matches = quad_to_grid(input, grid_size2, 0);
          //for (int m = 0; m < matches.size(); ++m) {
            //FloatType fraction = matches[m].fraction;
            //int index = matches[m].out;
            //int jj = index / grid_size_[2];
            //int ii = index % grid_size_[2];
            //for (int k = 0; k < shoebox_size_[0]; ++k) {
              //if (mask(k, j, i)) {
                //FloatType ivalue = image(k, j, i) * fraction;
                //FloatType bvalue = bkgrd(k, j, i) * fraction;
                //for (int kk = 0; kk < grid_size_[0]; ++kk) {
                  //FloatType zf = zfraction_(k, kk);
                  //profile_(kk, jj, ii) += ivalue * zf;
                  //background_(kk, jj, ii) += bvalue * zf;
                //}
              //}
            //}
          //}
        //}
      //}
    //}

    /**
     * Call the transform with the shoebox
     */
    //void call(const Shoebox<> &shoebox) {
      //af::versa< bool, af::c_grid<3> > mask(shoebox.mask.accessor());
      //af::ref< bool, af::c_grid<3> > mask_ref = mask.ref();
      //af::const_ref< int, af::c_grid<3> > temp_ref = shoebox.mask.const_ref();
      //for (std::size_t i = 0; i < mask_ref.size(); ++i) {
        //mask_ref[i] = (temp_ref[i] & Valid && temp_ref[i] & Foreground);
      //}
      //call(shoebox.data.const_ref(), shoebox.background.const_ref(), mask_ref);
    //}

    /**
     * Get a grid coordinate from an image coordinate
     * @param j The y index
     * @param i The x index
     * @returns The grid (c1, c2) index
     */
    //vec2<double> gc(std::size_t j, std::size_t i) const {
      //vec3<double> ds = s1_map_(y0_ + j, x0_ + i) - s1_;
      //return vec2<double>(grid_cent_[2] + (e1_ * ds) / step_size_[2],
                          //grid_cent_[1] + (e2_ * ds) / step_size_[1]);
    //}

    //int x0_, y0_;
    //int3 shoebox_size_;
    //int3 grid_size_;
    //double3 step_size_;
    //double3 grid_cent_;
    //vec3<double> s1_, e1_, e2_;
    //af::versa< vec3<double>, af::c_grid<2> > s1_map_arr_;
    //af::const_ref< vec3<double>, af::c_grid<2> > s1_map_;
    //af::versa< FloatType, af::c_grid<3> > profile_;
    //af::versa< FloatType, af::c_grid<3> > background_;
    //af::versa< FloatType, af::c_grid<2> > zfraction_arr_;
    //af::const_ref< FloatType, af::c_grid<2> > zfraction_;
  //};

}}}} // namespace dials::algorithms::reflection_basis::transform

#endif /* DIALS_ALGORITHMS_REFLEXION_BASIS_MAP_PIXELS_H */
