#ifndef DIALS_ALGORITHMS_REFLEXION_BASIS_REBIN_PIXELS_H
#define DIALS_ALGORITHMS_REFLEXION_BASIS_REBIN_PIXELS_H

namespace dials { namespace algorithms { namespace reflexion_basis {

  void rebin_pixels(const flex_double &output, const flex_double &input,
                    const flex_vec2_double &inputxy) {

    // Check the sizes
    DIALS_ASSERT(output.accessor().all().size() == 2);
    DIALS_ASSERT(input.accessor().all().size() == 2);
    DIALS_ASSERT(inputxy.accessor().all().size() == 2);
    DIALS_ASSERT(inputxy.accessor().all()[0] == input.accessor().all()[0] + 1);
    DIALS_ASSERT(inputxy.accessor().all()[1] == input.accessor().all()[1] + 1);

    // Get the input and output sizes
    std::size_t output_height = output.accessor().all()[0];
    std::size_t output_width = output.accessor().all()[1];
    std::size_t input_height = input.accessor().all()[0];
    std::size_t input_width = input.accessor().all()[1];

    // Loop through all the input pixels
    for (std::size_t j = 0; j < input_height; ++j) {
      for (std::size_t i = 0; i < input_width; ++i) {

        // Get the x, y coords of the target point
        ixy00 = inputxy(j, i);
        ixy01 = inputxy(j, i+1);
        ixy10 = inputxy(j+1, i);
        ixy11 = inputxy(j+1, i+1);

        // Create the target polygon and calculate its area
        vert4 target(ixy00, ixy01, ixy11, ixy10);
        double target_area = polygon_area(target);
        double value = input(j, i);

        // Get the range of new grid points
        double4 ix(ixy00[0], ixy01[0], ixy10[0], ixy11[0]);
        double4 iy(ixy00[1], ixy01[1], ixy10[1], ixy11[1]);
        int ox0 = (int)floor(min(ix)), ox1 = (int)ceil(max(ix));
        int oy0 = (int)floor(min(iy)), oy1 = (int)ceil(max(iy));

        // Cap the coordinates within the the output grid
        if (ox0 < 0) ox0 = 0;
        if (oy0 < 0) oy0 = 0;
        if (ox1 >= output_width) ox1 = output_width - 1;
        if (oy1 >= output_height) oy1 = output_height - 1;

        // Loop over all the pixels within the pixel range
        for (std::size_t jj = oy0; jj < oy1; ++jj) {
          for (std::size_t ii = ox0; ii < ox1; ++ii) {

            // Create the subject polygon
            vert4 subject(vec2<double>(ii, jj),
                          vec2<double>(ii+1, jj),
                          vec2<double>(ii+1,jj+1),
                          vec2<double>(ii,jj+1));

            // clip the polygon with the target polygon and calculate the
            // fraction of the area of the clipped polygon against the target.
            // Then redistribute the values from the target grid to the subject.
            vert8 result = quad_with_convex_quad(subject, target);
            double result_area = polygon_area(result);
            double fraction = result_area / target_area;
            output(jj, ii) += fraction * value;
          }
        }
      }
    }
  }

}}} // dials::algorithms::reflexion_basis

#endif /* DIALS_ALGORITHMS_REFLEXION_BASIS_REBIN_PIXELS_H */
