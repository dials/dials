
#ifndef DIALS_GEOMETRY_ALGORITHM_XDS_TRANSFORM_H
#define DIALS_GEOMETRY_ALGORITHM_XDS_TRANSFORM_H

#include <cmath>
#include <algorithm>
#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include <scitbx/array_family/flex_types.h>
#include "../transform/from_detector_to_xds.h"
#include "../detector_coordinate_system.h"
#include "../../equipment/detector.h"
#include "../../equipment/beam.h"
#include "../../equipment/goniometer.h"
#include "xds_transform_grid.h"

namespace dials { namespace geometry { namespace algorithm {

/** An array of vec3 <double> elements */
typedef scitbx::af::flex<scitbx::vec3 <double> >::type flex_vec3_double;

/** 
 * A class to calculate calculate the fraction of counts contributed by each 
 * data frame, j, around the reflection to each grid point, v3 in the profile 
 * frame.
 */
class xds_transform_e3_fraction {

public:

    /** Default constructor */
    xds_transform_e3_fraction() {}

    /** 
     * Initialise the fraction calculation class 
     * @param roi_size_z The image roi size in the z direction
     * @param grid_origin_e3 The XDS grid origin in the e3 direction
     * @param step_size_e3 The XDS grid step size in the e3 direction
     * @param starting_angle The rotation starting angle
     * @param oscillation_range The rotation oscillation range
     * @param sigma_mosaicity The standard deviation of the mosaicity
     */
    xds_transform_e3_fraction(int roi_size_z,
                              int grid_origin_e3,
                              double step_size_e3,
                              double starting_angle,
                              double oscillation_range,
                              double sigma_mosaicity) 
        : _roi_size_z(roi_size_z),
          _grid_origin_e3(grid_origin_e3),
          _step_size_e3(step_size_e3),
          _starting_angle(starting_angle),
          _oscillation_range(oscillation_range),
          _sigma_mosaicity(sigma_mosaicity) {}

    scitbx::af::flex_double calculate(double frame, double phi, double zeta);

private:

    int _roi_size_z;
    int _grid_origin_e3;
    double _step_size_e3;
    double _starting_angle;
    double _oscillation_range;
    double _sigma_mosaicity;
};

/**
 * Calculate the fraction of counts contributed by each data frame, j,
 * around the reflection to each grid point, v3 in the profile frame.
 *       
 * First we find and integrate over the set of phi angles covered by each 
 * data frame around the reflection to get Ij. Then we find the set of phi
 * angles covered by each grid coordinate. We then integrate over the 
 *
 * intersection of the phi angles covered by each data frame and each 
 * grid point to get Iv3j. The fraction of the counts is then calculated as 
 * Iv3j / Ij.
 *       
 * Further details of the method can be found in Kabsch 2010.
 *       
 * @param frame The z coordinate of the reflection (i.e. frame number)
 * @param phi The rotation angle of the reflection
 * @param zeta The lorentz correction factor
 *  
 * @returns An array containing the count fractions. The fraction of counts
 *          given by frame j to grid coordinate v3 can be found in the array
 *          by fv3j[v3-v30, j-j0]     
 *
 * @throws std::runtime_error if the supplied values are bad  
 */
scitbx::af::flex_double xds_transform_e3_fraction::calculate(
        double frame, double phi, double zeta)
{
    // Check the value of zeta
    if (zeta == 0.0) 
        throw std::runtime_error("zeta == 0.0"); 
        
    // Create an array to contain the intensity fractions
    int array_size = (2 * _grid_origin_e3 + 1) * (2 * _roi_size_z + 1);
    scitbx::af::flex_double fraction(array_size);

    // The range of data frames and grid points to iterate over
    int j0 = (int)frame - _roi_size_z;
    int j1 = (int)frame + _roi_size_z;
    int v30 = - _grid_origin_e3;
    int v31 = + _grid_origin_e3;
    
    // A constant used in the solution to the integrals below. 
    double sigr2 = 1.0 / (std::sqrt(2.0) * (_sigma_mosaicity / 
                                            std::abs(zeta)));

    // Loop over all j data frames in the region around the reflection
    for (int i = 0, j = j0; j <= j1; ++j) {

        // The data frame j covers the range of phi such that
        // rj = {phi':phi0 + (j-1)dphi <= phi' >= phi0 + jdpi}
        // Therefore the range of phi for j is given as follows.
        double bj = _starting_angle + j * _oscillation_range;
        double aj = bj - _oscillation_range;
        
        // Calculate the integral over rj (leaving out scaling factors):
        // I[exp(-(phi' - phi)^2 / (2 sigma^2)]
        double integral_j = (erf((bj - phi) * sigr2) - 
                             erf((aj - phi) * sigr2));
        
        // If integral is zero then set fractions to 0.0
        if (integral_j == 0.0) {
            for (int v3 = v30; v3 <= v31; ++v3, ++i) {
                fraction[i] = 0.0;
            }
        } else {
            double integral_j_r = 1.0 / integral_j;
            
            // Loop over all v3 in the profile grid
            for (int v3 = v30; v3 <= v31; ++v3) {

                // The grid coordinate v3 cover the range phi such that
                // rv3 = {phi':(v3 - 0.5)d3 <= (phi' - phi)zeta <= (v3 + 0.5)d3}
                // Therefore the range of phi for v3 is given as follows.
                double bv3 = ((v3 - 0.5) * _step_size_e3) / zeta + phi;
                double av3 = ((v3 + 0.5) * _step_size_e3) / zeta + phi;
                if (av3 > bv3) std::swap(av3, bv3);
                
                // We need to integrate over the intersection of sets rv3 and rj
                double av3j = std::max(av3, aj);
                double bv3j = std::min(bv3, bj);

                // If there is no intersection then set the fraction of the 
                // counts contributed by fata frame j to grid coordinate v3 to 
                // zero, otherwise calculate it as the ratio of the integral
                // over the intersection or rv3 and rj to the integral over rj
                if (av3j >= bv3j) {
                    fraction[i] = 0;
                } else {
                    fraction[i] = (erf((bv3j - phi) * sigr2) - 
                                   erf((av3j - phi) * sigr2)) * integral_j_r;
                }
                
                // Increment array index
                i++;
            }
        }
    }
    
    // Return the intensity fractions
    return fraction;
}

/**
 * A class to calculate the beam vector at each sub-divided pixel coordinate
 * that will be used during the transformation of the reflections from the
 * detector to the xds coordinate frame. This class is used during pre-
 * processing since no knowledge of the specific reflections are needed in order
 * to calculate the beam vectors. The beam vectors are then used along with
 * reflection specific stuff to calculate the xds coordinate for each pixel.
 */
class xds_transform_detector_beam_vectors {

public:

    /** The default constructor */
    xds_transform_detector_beam_vectors() {}    
    
    /**
     * Initialise the class.
     * @param detector The detector struct
     * @param wavelength The wavelength of the radiation
     * @param n_div The number of pixel sub-divisions
     */
    xds_transform_detector_beam_vectors(equipment::detector detector,
                                        double wavelength,
                                        int n_div)
        : _size(detector.get_size()),
          _origin(detector.get_origin()),
          _x_axis(detector.get_x_axis().normalize() * 
                detector.get_pixel_size()[0]),
          _y_axis(detector.get_y_axis().normalize() * 
                detector.get_pixel_size()[1]),
          _normal(detector.get_normal()),
          _distance(detector.get_distance()),
          _wavelength_r(1.0 / wavelength),
          _n_div(n_div) {}

    /**
     * Calculate the beam vector at every pixel on the detector, sub-divided
     * into (n_div * n_div) equal areas. This is done to remove a certain
     * amount of processing from being done per reflection and ensuring it
     * is only done before the reflections are procesed.
     * @returns An array of beam vectors 
     */
    flex_vec3_double calculate() {
    
        // Calculate the image size
        int x_size = _size[0] * _n_div;
        int y_size = _size[1] * _n_div;
        int image_size = x_size * y_size;
        double n_div_r = 1.0 / (double)_n_div;
        
        // Create the necessary arrays
        flex_vec3_double detector_s1x(x_size);
        flex_vec3_double detector_s1y(y_size);
        flex_vec3_double detector_s1(image_size);
        
        // Calculate the x and y components of the detector s1 vectors
        for (std::size_t i = 0; i < x_size; ++i) {
            detector_s1x[i] = (i * n_div_r - _origin[0]) * _x_axis;
        }
        for (std::size_t j = 0; j < y_size; ++j) {
            detector_s1y[j] = (j * n_div_r - _origin[1]) * _y_axis;
        }

        // Calculate the s1 vector for each sub-division of the detector        
        for (std::size_t k = 0, j = 0; j < y_size; ++j) {
            for (std::size_t i = 0; i < x_size; ++i, ++k) {
                detector_s1[k] = (
                    detector_s1x[i] + 
                    detector_s1y[j] + 
                    _distance * _normal).normalize() * _wavelength_r;
            }
        }
        
        // Return the s1 vector
        return detector_s1;
    }

private:
    
    equipment::detector detector;
    scitbx::vec2 <int> _size;
    scitbx::vec2 <double> _origin;
    scitbx::vec3 <double> _x_axis;
    scitbx::vec3 <double> _y_axis;
    scitbx::vec3 <double> _normal;
    double _distance;
    double _wavelength_r;
    int _n_div;
};

/**
 *
 */
class xds_transform {

public:

    /** The default constructor */
    xds_transform() {}

    /**
     * Initialise the transform.
     * @param grid The transform grid container
     * @param image The raw image volume
     * @param image_size The size of the image volume
     * @param detector The detector struct
     * @param beam The beam struct
     * @param gonio The goniometer struct
     * @param roi_size The region of interest to select (default (4, 4, 1)
     * @param n_div The number of pixel sub divisions to use (default 5)
     */
    xds_transform(xds_transform_grid grid,
                  scitbx::af::flex_int image,
                  scitbx::vec3 <int> image_size,
                  equipment::detector detector,
                  equipment::beam beam,
                  equipment::goniometer gonio,
                  scitbx::vec3 <int> roi_size = scitbx::vec3 <int> (4, 4, 1),
                  int n_div = 5)
        : 
          _grid(grid),
          _image(image),
          _image_size(image_size),
          _roi_size(roi_size),
          _grid_origin(grid.get_origin()),
          _step_size(grid.get_step_size()),
          _grid_size(grid.get_size()),
          _starting_frame(gonio.get_starting_frame()),
          _m2(gonio.get_rotation_axis()),
          _s0(beam.get_direction()),
          _n_div(n_div)
    {
        // Create an object to calculate the detector beam vectors. Then 
        // Calculate and save the detector beam vectors.
        xds_transform_detector_beam_vectors beam_vectors(
                detector, beam.get_wavelength(), n_div);
        _detector_s1 = beam_vectors.calculate();
        
        // Initialise an object to calculate the e3 fraction vector
        _e3_fraction = xds_transform_e3_fraction(
                roi_size[2], grid.get_origin()[2], grid.get_step_size()[2], 
                gonio.get_starting_angle(), gonio.get_oscillation_range(), 
                grid.get_sigma_mosaicity());
    }

    void calculate(scitbx::vec3 <double> xyz, 
                   scitbx::vec3 <double> s1,
                   double phi);

private:

    xds_transform_e3_fraction _e3_fraction;
    xds_transform_grid _grid;
    flex_vec3_double _detector_s1;
    scitbx::af::flex_int _image;
    scitbx::vec3 <int> _image_size;
    scitbx::vec3 <int> _roi_size;
    scitbx::vec3 <int> _grid_size;
    scitbx::vec3 <int> _grid_origin;
    scitbx::vec3 <double> _step_size;
    scitbx::vec3 <double> _s0;
    scitbx::vec3 <double> _m2;
    double _starting_frame;
    int _n_div;
};

/**
 * @param xyz The coordinate of the reflection in the detector image volume
 * @param s1 The beam vector of the reflection
 * @param phi The rotation angle of the reflection
 */
void xds_transform::calculate(scitbx::vec3 <double> xyz, 
                              scitbx::vec3 <double> s1, 
                              double phi)
{
    // Constant for scaling values
    static const double r2d = 1.0 / scitbx::constants::pi_180;

    int x_stride = _image_size[0] * _n_div;
    int y_stride = _image_size[1] * _n_div;

    // Get the grid data array
    scitbx::af::flex_double grid = _grid.get_data();

    // Calculate the x, y, z ranges to iterate over
    int x0 = ((int)xyz[0] - _roi_size[0]) * _n_div;
    int x1 = ((int)xyz[0] + _roi_size[0]) * _n_div;
    int y0 = ((int)xyz[1] - _roi_size[1]) * _n_div;
    int y1 = ((int)xyz[1] + _roi_size[1]) * _n_div;
    int z0 = ((int)xyz[2] - _roi_size[2] - _starting_frame);
    int z1 = ((int)xyz[2] + _roi_size[2] - _starting_frame); 
    
    // Calculate 1 / n_div and 1 / (n_div*n_div) for convenience
    double n_div_r = 1.0 / _n_div;
    double div_fraction = n_div_r * n_div_r;

    // Calculate the reflection coordinate system e1 and e2 axes, and zeta, the
    // loretz correction (used to calculate component on e3 axis
    scitbx::vec3 <double> e1 = s1.cross(_s0).normalize();
    scitbx::vec3 <double> e2 = s1.cross(e1).normalize();
    double zeta = _m2 * e1;
    double s1_length = s1.length();
    e1 = e1 * r2d / s1_length;
    e2 = e2 * r2d / s1_length;

    // Calculate e1.s1 and e2.s1 here. The e1 and e2 coordinates are calculated
    // as e1.(s' - s1) and e2.(s' - s1). Putting this here means we add an extra
    // 6 multiplications but remove 2 * nx * ny subtractions. Only really saves
    // a small amount of time.
    double c11 = e1 * s1;
    double c21 = e2 * s1;

    // Calculate the fraction of counts contributed by each data frame, j, 
    // around the reflection to each grid point, v3 in the profile frame. Hold
    // these fractions in a 2d array.
    scitbx::af::flex_double fraction = _e3_fraction.calculate(xyz[2], phi, zeta); 
    
    // Loop through all the pixels (and their sub-divisions). Calculate the
    // coordinate of each pixel in the XDS coordinate frame e1 and e2 axes.
    // Find the grid point in which the calculate point is contained and then
    // add the counts for that pixel to the grid. See Kabsch 2010
    for (int yy = y0; yy < y1; ++yy) {
        for (int xx = x0; xx < x1; ++xx) {
            double c1 = e1 * _detector_s1[xx + yy * x_stride] - c11;
            double c2 = e2 * _detector_s1[xx + yy * x_stride] - c21;
            int gi = _grid_origin[0] + c1 / _step_size[0];
            int gj = _grid_origin[1] + c2 / _step_size[1];
            if (gi < 0 || gi >= 9 || gj < 0 || gj >= 9) {
                continue;
            }
            int x = xx * n_div_r;
            int y = yy * n_div_r;
            for (int z = z0; z <= z1; ++z) {
                
                int image_index = z * _image_size[1] * _image_size[0] +
                                  y * _image_size[0] +
                                  x;
                int value = _image[image_index];               
                
                for (int gk = 0; gk < _grid_size[2]; ++gk) {
                    int grid_index = gk * _grid_size[1] * _grid_size[0] +
                                     gj * _grid_size[0] +
                                     gi;
                    int fraction_index = (z - z0) * _grid_size[2] + gk;
                    grid[grid_index] += value * fraction[fraction_index] * div_fraction;
                }                        
            }
        }
    }
}

}}}

#endif // DIALS_GEOMETRY_ALGORITHM_XDS_TRANSFORM_H
