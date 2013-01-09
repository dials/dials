
#ifndef DIALS_GEOMETRY_ALGORITHM_XDS_TRANSFORM_H
#define DIALS_GEOMETRY_ALGORITHM_XDS_TRANSFORM_H

#include <cmath>
#include <algorithm>
#include <scitbx/vec3.h>
#include <scitbx/array_family/flex_types.h>
#include "../transform/from_detector_to_xds.h"

namespace dials { namespace geometry { namespace algorithm {

class xds_transform_e3_fraction {

public:

    xds_transform_e3_fraction() {}

    xds_transform_e3_fraction(int roi_size_z,
                              int grid_size_e3,
                              double step_size_e3,
                              double starting_angle,
                              double oscillation_range,
                              double sigma_mosaicity) 
        : _roi_size_z(roi_size_z),
          _grid_size_e3(grid_size_e3),
          _step_size_e3(step_size_e3),
          _starting_angle(starting_angle),
          _oscillation_range(oscillation_range),
          _sigma_mosaicity(sigma_mosaicity) {}

    scitbx::af::flex_double calculate(double frame, double phi, double zeta);

private:

    int _roi_size_z;
    int _grid_size_e3;
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
 * @param rcs The reflection coordinate system
 * @param frame The z coordinate of the reflection (i.e. frame number)
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
    int array_size = (2 * _grid_size_e3 + 1) * (2 * _roi_size_z + 1);
    scitbx::af::flex_double fraction(array_size);

    // The range of data frames and grid points to iterate over
    int j0 = (int)frame - _roi_size_z;
    int j1 = (int)frame + _roi_size_z;
    int v30 = - _grid_size_e3;
    int v31 = + _grid_size_e3;
    
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


class xds_transform {

public:

    xds_transform() {}

    xds_transform(scitbx::af::flex_int image,
                  scitbx::vec3 <int> image_size,
                  scitbx::vec3 <int> roi_size,
                  scitbx::vec3 <int> grid_size,
                  scitbx::vec3 <double> step_size,
                  double starting_frame,
                  double starting_angle,
                  double oscillation_range,
                  double sigma_mosaicity) 
        : _e3_fraction(roi_size[2], 
                       grid_size[2], 
                       step_size[2], 
                       starting_angle, 
                       oscillation_range, 
                       sigma_mosaicity),
          _image(image),
          _image_size(image_size),
          _roi_size(roi_size),
          _grid_size(grid_size),
          _step_size(step_size),
          _starting_frame(starting_frame),
          _starting_angle(starting_angle),
          _oscillation_range(oscillation_range) {}

    scitbx::af::flex_double calculate(double frame, double phi, double zeta);

    scitbx::af::flex_double calculate2(transform::from_detector_to_xds dcs_to_xcs,
                    scitbx::vec3 <double> xyz, double phi, double zeta);

    scitbx::af::flex_double calculate3(transform::from_detector_to_xds dcs_to_xcs,
                    scitbx::vec3 <double> xyz, double phi, double zeta);

private:

    xds_transform_e3_fraction _e3_fraction;
    scitbx::af::flex_int _image;
    scitbx::vec3 <int> _image_size;
    scitbx::vec3 <int> _roi_size;
    scitbx::vec3 <int> _grid_size;
    scitbx::vec3 <double> _step_size;
    double _starting_frame;
    double _starting_angle;
    double _oscillation_range;
};

scitbx::af::flex_double xds_transform::calculate(double frame, double phi, double zeta)
{
    return _e3_fraction.calculate(frame, phi, zeta);
}  


scitbx::af::flex_double xds_transform::calculate3(
                      transform::from_detector_to_xds dcs_to_xcs,
                      scitbx::vec3 <double> xyz, double phi, double zeta)
{
    scitbx::vec3 <int> grid_dim = 2 * _grid_size + 1;
    int grid_array_size = grid_dim[0] * grid_dim[1] * grid_dim[2];  

    scitbx::af::flex_double grid(grid_array_size);

    int x0 = (int)xyz[0] - _roi_size[0];
    int x1 = (int)xyz[0] + _roi_size[0];
    int y0 = (int)xyz[1] - _roi_size[1];
    int y1 = (int)xyz[1] + _roi_size[1];
    int z0 = (int)xyz[2] - _roi_size[2] - _starting_frame;
    int z1 = (int)xyz[2] + _roi_size[2] - _starting_frame; 
    
    int n_div = 5;
    double n_div_r = 1.0 / n_div;
    double div_fraction = n_div_r * n_div_r;
    
    
    scitbx::af::flex_double fraction = _e3_fraction.calculate(xyz[2], phi, zeta); 
    
    for (int y = y0; y <= y1; ++y) {
        for (int x = x0; x <= x1; ++x) {
            for (int z = z0; z <= z1; ++z) {
                
                int image_index = z * _image_size[1] * _image_size[0] +
                                  y * _image_size[0] +
                                  x;
                int value = _image[image_index];                
                
                for (int yy = 0; yy < n_div; ++yy) {
                    for (int xx = 0; xx < n_div; ++xx) {
                        double xxx = x + xx * n_div_r;
                        double yyy = y + yy * n_div_r;

                        double phi_dash = _starting_angle + (z - _starting_frame) 
                            * _oscillation_range;

                        scitbx::vec3 <double> c = dcs_to_xcs.apply(
                            scitbx::vec2 <double> (xxx, yyy), phi_dash);                     
                        
                        int gi = _grid_size[0] + c[0] / _step_size[0];
                        int gj = _grid_size[1] + c[1] / _step_size[1];
                        
                        for (int gk = 0; gk < grid_dim[2]; ++gk) {
                                                        
                            if (0 <= gi && gi < 9 &&
                                0 <= gj && gj < 9 &&
                                0 <= gk && gk < 9) {
                                int grid_index = gk * grid_dim[1] * grid_dim[0] +
                                                 gj * grid_dim[0] +
                                                 gi;
                                int fraction_index = (z - z0) * grid_dim[2] + gk;
                                grid[grid_index] += value * fraction[fraction_index] * div_fraction;
                            }
                        }                        
                    }
                }
                
            }
        }
    }
    
    return grid;
}

scitbx::af::flex_double xds_transform::calculate2(
                      transform::from_detector_to_xds dcs_to_xcs,
                      scitbx::vec3 <double> xyz, double phi, double zeta)
{
    // Calculate the range of frames to use
    int x0 = (int)xyz[0] - _roi_size[0];
    int x1 = (int)xyz[0] + _roi_size[0];
    int y0 = (int)xyz[1] - _roi_size[1];
    int y1 = (int)xyz[1] + _roi_size[1];
    int z0 = (int)xyz[2] - _roi_size[2];
    int z1 = (int)xyz[2] + _roi_size[2]; 

    // Calculate the e3 fraction
    scitbx::af::flex_double fraction = _e3_fraction.calculate(xyz[2], phi, zeta);    
    
    int grid_size = (2 * _grid_size[0] + 1) *
                    (2 * _grid_size[1] + 1) *
                    (2 * _grid_size[2] + 1);
    scitbx::vec3 <int> grid_dim = 2 * _grid_size + 1;
    scitbx::af::flex_double grid(grid_size);
    
    int n_div = 5;
    double n_div_r = 1.0 / (double)n_div;
    double div_fraction_r = 1.0 / (double)(n_div * n_div);

    for (int j = y0; j <= y1; ++j) {
        for (int i = x0; i <= x1; ++i) {
        
            for (double jj = j; jj < j + 1; jj += n_div_r) {
                for (double ii = i; ii < i + 1; ii += n_div_r) {
              
                    for (int k = z0; k <= z1; ++k) {
                        
                        double phi_dash = _starting_angle + (k - _starting_frame) 
                            * _oscillation_range;
                
                        scitbx::vec3 <double> c = dcs_to_xcs.apply(
                            scitbx::vec2 <double> (ii, jj), phi_dash);  
              
                        std::cout << c[0] << " " << c[1] << " " << c[2] << std::endl;
                        
                        int iindex = k * (_image_size[0]*_image_size[1]) + j * _image_size[0] + i;
                        int value = _image[iindex];
                        
                        int gi = _grid_size[0] + (int)(c[0] / _step_size[0]);
                        int gj = _grid_size[1] + (int)(c[1] / _step_size[1]);
                        int gk = _grid_size[2] + (int)(c[2] / _step_size[2]);
                        
                        std::cout << gi << " " << gj << " " << gk << std::endl;
                        
                        if (0 <= gi && gi < 2 * _grid_size[1] + 1 &&
                            0 <= gj && gj < 2 * _grid_size[2] + 1 &&
                            0 <= gk && gk < 2 * _grid_size[2] + 1) {
                            int gindex = gk * (grid_dim[0] * grid_dim[1]) + gj * grid_dim[0] + gi;
                            grid[gindex] += value * fraction[k] * div_fraction_r;
                        }
                    }
                }
            }        
        }        
    }
    
    return grid;
}

}}}

#endif // DIALS_GEOMETRY_ALGORITHM_XDS_TRANSFORM_H
