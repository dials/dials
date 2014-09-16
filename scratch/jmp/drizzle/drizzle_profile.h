
#ifndef DIALS_ALGORITHMS_INTEGRATION_DRIZZLE_PROFILE_H
#define DIALS_ALGORITHMS_INTEGRATION_DRIZZLE_PROFILE_H

namespace dials { namespace algorithms {

  class ProfileDrizzler {
  public:

    class Transform {
    public:

      Transform(int4 bbox, const CoordinateSystem &cs)
        : bbox_(bbox),
          cs_(cs) {}

      vec2<double> operator()(vec2<double> xy) const {
        xy[0] += bbox_[0];
        xy[1] += bbox_[2];
        vec3<double> s1 = detector_.get_pixel_lab_coord(xy);
        vec2<double> e12 = cs_.from_beam_vector(s1);
      }

    private:

      int bbox_;
      const CoordinateSystem &cs_;
    };

    ProfileDrizzler(int3 grid_size, double fraction)
      : profile_(grid_size),
        profile_weight_(grid_size),
        fraction_(fraction) {

    }

    af::versa< double, af::c_grid<2> > profile() const {
      return profile_;
    }

    af::versa< double, af::c_grid<2> > weight() const {
      return profile_weight_;
    }

    void add(const af::const_ref< double, af::c_grid<2> > &image,
             const af::const_ref< double, af::c_grid<2> > &weight,
             const CoordinateSystem &cs) {
      drizzle(
          profile_.ref(),
          profile_weight_.ref(),
          image,
          weight,
          fraction_,
          make_transform(cs));
    }

    void add(const af::const_ref< double, af::c_grid<2> > &image,
             const CoordinateSystem &cs) {
      af::versa< double, af::c_grid<2> > weight(image.accessor(), 1.0);
      add(image, weight.const_ref(), cs);
    }

  private:

    Transform make_transform(const CoordinateSystem &cs) const {
      return Transform(cs);
    }

    af::versa< double, af::c_grid<2> > profile_;
    af::versa< double, af::c_grid<2> > profile_weight_;
    double fraction_;
  };

}}

#endif // DIALS_ALGORITHMS_INTEGRATION_DRIZZLE_PROFILE_H
