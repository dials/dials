
#ifndef DIALS_ALGORITHMS_IMAGE_DRIZZLE_DRIZZLE_H
#define DIALS_ALGORITHMS_IMAGE_DRIZZLE_DRIZZLE_H

#include <dials/error.h>

namespace dials { namespace algorithms {

  template <typename T>
  T min4(T a, T b, T c, T d) {
    return std::min(std::min(a,b),std::min(c,d));
  }

  template <typename T>
  T max4(T a, T b, T c, T d) {
    return std::max(std::max(a,b),std::max(c,d));
  }

  template <typename Transform>
  void drizzle(af::ref<double, af::c_grid<2> > dst,
               af::ref<double, af::c_grid<2> > dst_weight,
               const af::const_ref<double, af::c_grid<2> > &src,
               const af::const_ref<double, af::c_grid<2> > &src_weight,
               double fraction,
               Transform transform) {

    // Ensure the input is valid
    DIALS_ASSERT(fraction > 0 && fraction <= 1);
    DIALS_ASSERT(dst.accessor().all_eq(dst_weight.accessor()));
    DIALS_ASSERT(src.accessor().all_eq(src_weight.accessor()));

    // The offset into the pixel
    double d = fraction / 2;

    // The size of the output grid
    std::size_t height = dst.accessor()[0];
    std::size_t width = dst.accessor()[1];

    // Loop over the src pixels
    for (std::size_t j = 0; j < src.accessor()[0]; ++j) {
      for (std::size_t i = 0; i < src.accessor()[1]; ++i) {

        // The smaller internal pixel
        double i0 = i + 0.5 - d;
        double i1 = i + 0.5 + d;
        double j0 = j + 0.5 - d;
        double j1 = j + 0.5 + d;

        // Transform the four points
        vec2<double> p1 = transform(vec2<double>(i0, j0));
        vec2<double> p2 = transform(vec2<double>(i1, j0));
        vec2<double> p3 = transform(vec2<double>(i0, j1));
        vec2<double> p4 = transform(vec2<double>(i1, j1));
        vert4 a(p1, p2, p3, p4);

        // The range this pixel covers on the output grid,
        // clipped by the size of the output grid
        int x0 = min4(p1[0], p2[0], p3[0], p4[0]);
        int x1 = max4(p1[0], p2[0], p3[0], p4[0]);
        int y0 = min4(p1[1], p2[1], p3[1], p4[1]);
        int y1 = max4(p1[1], p2[1], p3[1], p4[1]);
        x0 = std::max(x0, 0);
        y0 = std::max(y0, 0);
        x1 = std::min(x1, (int)width);
        y1 = std::min(y1, (int)height);

        // Transformed area
        double total_area = simple_area(a);

        // Loop over possible intersections and compute fraction of intensity
        // contributed to each output grid point
        for (int y = y0; y < y1; ++y) {
          for (int x = x0; x < x1; ++x) {
            vert4 b(vec2<double>(x,y),
                    vec2<double>(x+1,y),
                    vec2<double>(x+1,y+1),
                    vec2<double>(x,y+1))
            double frac_area = simple_area(quad_with_convex_quad(a, b));
            double frac_weight = frac_area * src_weight(j,i);
            double weight = frac_weight + dst_weight(y,x);
            if (weight != 0) {
              double value = frac_weight*src(j,i) + dst(y,x)*dst_weight(y,x);
              dst(y,x) = num / weight;
            }
            dst_weight(y,x) = weight;
          }
        }
      }
    }
  }

}}

#endif // DIALS_ALGORITHMS_IMAGE_DRIZZLE_DRIZZLE_H
