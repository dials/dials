/*
 * map_pixels.h
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
#include <dials/algorithms/polygon/spatial_interpolation.h>
#include <dials/algorithms/reflection_basis/coordinate_system.h>
#include <dials/model/data/shoebox.h>
#include <dials/model/data/reflection.h>

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
  using dials::model::Reflection;
  using dials::model::Foreground;
  using dials::model::Valid;
  using dials::algorithms::polygon::spatial_interpolation::Match;
  using dials::algorithms::polygon::spatial_interpolation::quad_to_grid;

  class TransformSpec {
  public:
    TransformSpec(const Beam &beam, const Detector &detector,
                  const Goniometer &gonio, const Scan &scan,
                  double mosaicity, double n_sigma, std::size_t grid_size)
      : s0_(beam.get_s0()),
        m2_(gonio.get_rotation_axis().normalize()),
        image_size_(detector.get_image_size()[1], detector.get_image_size()[0]),
        grid_size_(2*grid_size+1, 2*grid_size+1, 2*grid_size+1),
        step_size_(mosaicity * n_sigma / grid_size,
                   beam.get_sigma_divergence() * n_sigma / grid_size,
                   beam.get_sigma_divergence() * n_sigma / grid_size),
        grid_centre_(grid_size + 0.5, grid_size + 0.5, grid_size + 0.5),
        s1_map_(beam_vector_map(detector, beam, true)),
        map_frames_(scan.get_oscillation()[0],
                    scan.get_oscillation()[1],
                    mosaicity, n_sigma, grid_size) {
      DIALS_ASSERT(image_size_.all_gt(0));
      DIALS_ASSERT(step_size_.all_gt(0));
      DIALS_ASSERT(grid_size_.all_gt(0));
    }

    vec3<double> m2() const {
      return m2_;
    }

    vec3<double> s0() const {
      return s0_;
    }

    int2 image_size() const {
      return image_size_;
    }

    int3 grid_size() const {
      return grid_size_;
    }

    double3 step_size() const {
      return step_size_;
    }

    double3 grid_centre() const {
      return grid_centre_;
    }

    af::versa< vec3<double>, af::c_grid<2> > s1_map() const {
      return s1_map_;
    }

    af::versa< double, af::c_grid<2> > map_frames(
        vec2<int> frames, double phi, double zeta) const {
      return map_frames_(frames, phi, zeta);
    }

  private:
    vec3<double> s0_;
    vec3<double> m2_;
    int2 image_size_;
    int3 grid_size_;
    double3 step_size_;
    double3 grid_centre_;
    af::versa< vec3<double>, af::c_grid<2> > s1_map_;
    MapFramesForward map_frames_;
  };


  class Forward2 {
  public:

    Forward2(const TransformSpec &spec,
            const vec3<double> &s1, double phi, int6 bbox,
            const af::const_ref< double, af::c_grid<3> > &image,
            const af::const_ref< bool, af::c_grid<3> > &mask) {
      init(spec, s1, phi, bbox);
      call(image, mask);
    }

    Forward2(const TransformSpec &spec,
            const vec3<double> &s1, double phi, int6 bbox,
            const af::const_ref< double, af::c_grid<3> > &image,
            const af::const_ref< double, af::c_grid<3> > &bkgrd,
            const af::const_ref< bool, af::c_grid<3> > &mask) {
      init(spec, s1, phi, bbox);
      call(image, bkgrd, mask);
    }

    Forward2(const TransformSpec &spec,
            const CoordinateSystem &cs, int6 bbox,
            const af::const_ref< double, af::c_grid<3> > &image,
            const af::const_ref< bool, af::c_grid<3> > &mask) {
      init(spec, cs, bbox);
      call(image, mask);
    }

    Forward2(const TransformSpec &spec,
            const CoordinateSystem &cs, int6 bbox,
            const af::const_ref< double, af::c_grid<3> > &image,
            const af::const_ref< double, af::c_grid<3> > &bkgrd,
            const af::const_ref< bool, af::c_grid<3> > &mask) {
      init(spec, cs, bbox);
      call(image, bkgrd, mask);
    }

    Forward2(const TransformSpec &spec,
            const Reflection &reflection) {
      init(spec, reflection);
      af::const_ref< int, af::c_grid<3> > shoebox_mask =
        reflection.get_shoebox_mask().const_ref();
      af::versa< bool, af::c_grid<3> > mask(shoebox_mask.accessor(), false);
      for (std::size_t i = 0; i < mask.size(); ++i) {
        mask[i] = shoebox_mask[i] & Valid && shoebox_mask[i] & Foreground;
      }
      call(reflection.get_shoebox().const_ref(),
           reflection.get_shoebox_background().const_ref(),
           mask.const_ref());
    }

    af::versa< double, af::c_grid<3> > profile() const {
      return profile_;
    }

    af::versa< double, af::c_grid<3> > background() const {
      return background_;
    }

    af::versa< double, af::c_grid<2> > zfraction() const {
      return zfraction_arr_;
    }

  private:

    void init(const TransformSpec &spec, const vec3<double> &s1,
              double phi, int6 bbox) {
      CoordinateSystem cs(spec.m2(), spec.s0(), s1, phi);
      init(spec, cs, bbox);
    }

    void init(const TransformSpec &spec, const Reflection &r) {
      vec3<double> s1 = r.get_beam_vector();
      double phi = r.get_rotation_angle();
      int6 bbox = r.get_bounding_box();
      CoordinateSystem cs(spec.m2(), spec.s0(), s1, phi);
      init(spec, cs, bbox);
    }

    void init(const TransformSpec &spec,
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

      // Calculate the fraction of intensity contributed from each data
      // frame to each grid coordinate
      vec2<int> zrange(bbox[4], bbox[5]);
      zfraction_arr_ = spec.map_frames(zrange, cs.phi(), cs.zeta());
      zfraction_ = zfraction_arr_.const_ref();
    }

    void call(const af::const_ref< double, af::c_grid<3> > &image,
              const af::const_ref< bool, af::c_grid<3> > &mask) {

      // Check the input
      DIALS_ASSERT(image.accessor().all_eq(shoebox_size_));
      DIALS_ASSERT(image.accessor().all_eq(mask.accessor()));

      // Initialise the profile arrays
      af::c_grid<3> accessor(grid_size_);
      profile_ = af::versa< double, af::c_grid<3> >(accessor, 0.0);

      // Loop through all the points in the shoebox. Calculate the polygon
      // formed by the pixel in the local coordinate system. Find the points
      // on the grid which intersect with the polygon and the fraction of the
      // pixel area shared with each grid point. For each intersection, loop
      // through the frames, mapping the fraction of the pixel value in each
      // frame to the grid point.
      int2 grid_size2(grid_size_[1], grid_size_[2]);
      for (std::size_t j = 0; j < shoebox_size_[1]; ++j) {
        for (std::size_t i = 0; i < shoebox_size_[2]; ++i) {
          vert4 input(gc(j, i), gc(j, i+1), gc(j+1, i+1), gc(j+1, i));
          af::shared<Match> matches = quad_to_grid(input, grid_size2, 0);
          for (int m = 0; m < matches.size(); ++m) {
            double fraction = matches[m].fraction;
            int index = matches[m].out;
            int ii = index % grid_size_[2];
            int jj = index / grid_size_[2];
            for (int k = 0; k < shoebox_size_[0]; ++k) {
              if (mask(k, j, i)) {
                double value = image(k, j, i) * fraction;
                for (int kk = 0; kk < grid_size_[0]; ++kk) {
                  profile_(kk, jj, ii) += value * zfraction_(k, kk);
                }
              }
            }
          }
        }
      }
    }

    void call(const af::const_ref< double, af::c_grid<3> > &image,
              const af::const_ref< double, af::c_grid<3> > &bkgrd,
              const af::const_ref< bool, af::c_grid<3> > &mask) {

      // Check the input
      DIALS_ASSERT(image.accessor().all_eq(shoebox_size_));
      DIALS_ASSERT(image.accessor().all_eq(mask.accessor()));
      DIALS_ASSERT(image.accessor().all_eq(bkgrd.accessor()));

      // Initialise the profile arrays
      af::c_grid<3> accessor(grid_size_);
      profile_ = af::versa< double, af::c_grid<3> >(accessor, 0.0);
      background_ = af::versa< double, af::c_grid<3> >(accessor, 0.0);

      // Loop through all the points in the shoebox. Calculate the polygon
      // formed by the pixel in the local coordinate system. Find the points
      // on the grid which intersect with the polygon and the fraction of the
      // pixel area shared with each grid point. For each intersection, loop
      // through the frames, mapping the fraction of the pixel value in each
      // frame to the grid point.
      int2 grid_size2(grid_size_[1], grid_size_[2]);
      for (std::size_t j = 0; j < shoebox_size_[1]; ++j) {
        for (std::size_t i = 0; i < shoebox_size_[2]; ++i) {
          vert4 input(gc(j, i), gc(j, i+1), gc(j+1, i+1), gc(j+1, i));
          af::shared<Match> matches = quad_to_grid(input, grid_size2, 0);
          for (int m = 0; m < matches.size(); ++m) {
            double fraction = matches[m].fraction;
            int index = matches[m].out;
            int jj = index / grid_size_[2];
            int ii = index % grid_size_[2];
            for (int k = 0; k < shoebox_size_[0]; ++k) {
              if (mask(k, j, i)) {
                double ivalue = image(k, j, i) * fraction;
                double bvalue = image(k, j, i) * fraction;
                for (int kk = 0; kk < grid_size_[0]; ++kk) {
                  double zf = zfraction_(k, kk);
                  profile_(kk, jj, ii) += ivalue * zf;
                  background_(kk, jj, ii) += bvalue * zf;
                }
              }
            }
          }
        }
      }
    }

    vec2<double> gc(std::size_t j, std::size_t i) const {
      vec3<double> ds = s1_map_(y0_ + j, x0_ + i) - s1_;
      return vec2<double>(grid_cent_[2] + (e1_ * ds) / step_size_[2],
                          grid_cent_[1] + (e2_ * ds) / step_size_[1]);
    }

    int x0_, y0_;
    int3 shoebox_size_;
    int3 grid_size_;
    double3 step_size_;
    double3 grid_cent_;
    vec3<double> s1_, e1_, e2_;
    af::versa< vec3<double>, af::c_grid<2> > s1_map_arr_;
    af::const_ref< vec3<double>, af::c_grid<2> > s1_map_;
    af::versa< double, af::c_grid<3> > profile_;
    af::versa< double, af::c_grid<3> > background_;
    af::versa< double, af::c_grid<2> > zfraction_arr_;
    af::const_ref< double, af::c_grid<2> > zfraction_;
  };




  class ForwardBatch {
  public:

    typedef af::versa< double, af::c_grid<3> > versa_double3;
    typedef af::versa< double, af::c_grid<4> > versa_double4;

    ForwardBatch(const TransformSpec &spec,
                 const af::const_ref<Reflection> &rlist)
      : size_(spec.grid_size()),
        num_(rlist.size()),
        profile_(init_profile_grid(), 0.0),
        background_(init_profile_grid(), 0.0) {
      std::size_t offset = 0;
      for (std::size_t i = 0; i < rlist.size(); ++i) {
        Forward2 transform(spec, rlist[i]);
        versa_double3 p = transform.profile();
        versa_double3 b = transform.background();
        std::copy(p.begin(), p.end(), profile_.begin() + offset);
        std::copy(b.begin(), b.end(), background_.begin() + offset);
        offset += p.size();
      }
    }

    versa_double4 profile() const {
      return profile_;
    }

    versa_double4 background() const {
      return background_;
    }

    versa_double3 profile(std::size_t index) const {
      DIALS_ASSERT(index < num_);
      versa_double3 p(af::c_grid<3>(size_[0], size_[1], size_[2]),
        af::init_functor_null<double>());
      std::size_t beg = index * p.size(), end = beg + p.size();
      std::copy(profile_.begin() + beg,  profile_.begin() + end, p.begin());
      return p;
    }

    versa_double3 background(std::size_t index) const {
      DIALS_ASSERT(index < num_);
      versa_double3 b(af::c_grid<3>(size_[0], size_[1], size_[2]),
        af::init_functor_null<double>());
      std::size_t beg = index * b.size(), end = beg + b.size();
      std::copy(background_.begin() + beg, background_.begin() + end, b.begin());
      return b;
    }

  private:

    af::c_grid<4> init_profile_grid() const {
      af::c_grid<4> grid;
      grid[0] = num_;
      grid[1] = size_[0];
      grid[2] = size_[1];
      grid[3] = size_[2];
      return grid;
    }

    int3 size_;
    std::size_t num_;
    versa_double4 profile_;
    versa_double4 background_;
  };







  inline
  void forward_batch(const TransformSpec &spec, af::ref<Reflection> rlist) {
    for (std::size_t i = 0; i < rlist.size(); ++i) {
      if (rlist[i].is_valid()) {
        Forward2 transform(spec, rlist[i]);
        rlist[i].set_transformed_shoebox(transform.profile());
  //      rlist[i].set_transformed_shoebox_background(transform.background());
      }
    }
  }

}}}} // namespace dials::algorithms::reflection_basis::transform

#endif /* DIALS_ALGORITHMS_REFLEXION_BASIS_MAP_PIXELS_H */
