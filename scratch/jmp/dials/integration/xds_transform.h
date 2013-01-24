
#ifndef DIALS_INTEGRATION_XDS_TRANSFORM_H
#define DIALS_INTEGRATION_XDS_TRANSFORM_H

#include <cmath>
#include <algorithm>
#include <exception>
#include <boost/math/special_functions/erf.hpp>
#include <scitbx/constants.h>
#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include <scitbx/array_family/flex_types.h>
#include "../error.h"
#include "../equipment/detector.h"
#include "../equipment/beam.h"
#include "../equipment/goniometer.h"
#include "../array_family/array_types.h"
#include "xds_transform_grid.h"

namespace dials { namespace integration {

/** 
 * A class to calculate calculate the fraction of counts contributed by each 
 * data frame, j, around the reflection to each grid point, v3 in the profile 
 * frame.
 */
class XdsTransformE3Fraction {

public:

    /** Default constructor */
    XdsTransformE3Fraction() {}

    /** 
     * Initialise the fraction calculation class 
     * @param goniometer The goniometer parameters
     * @param grid_origin_e3 The XDS grid origin in the e3 direction
     * @param step_size_e3 The XDS grid step size in the e3 direction
     * @param sigma_mosaicity The standard deviation of the mosaicity
     */
    XdsTransformE3Fraction(const equipment::Goniometer &goniometer,
                           const XdsTransformGrid &grid)
        : starting_angle_(goniometer.get_starting_angle()),
          oscillation_range_(goniometer.get_oscillation_range()),
          starting_frame_(goniometer.get_starting_frame()),
          grid_origin_e3_(grid.get_origin()[0]),
          step_size_e3_(grid.get_step_size()[0]),
          sigma_mosaicity_(grid.get_sigma_mosaicity()) {}

    scitbx::af::flex_double calculate(scitbx::vec2 <int> roi_z,
                                      double phi, double zeta);

private:

    double starting_angle_;
    double oscillation_range_;
    int starting_frame_;
    int grid_origin_e3_;
    double step_size_e3_;
    double sigma_mosaicity_;
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
scitbx::af::flex_double XdsTransformE3Fraction::calculate(
        scitbx::vec2 <int> roi_z, double phi, double zeta)
{
    // Check the value of zeta
    DIALS_ASSERT(zeta != 0.0);
        
    // Create an array to contain the intensity fractions
    int array_size = (2 * grid_origin_e3_ + 1) * (roi_z[1] + 1 - roi_z[0]);
    scitbx::af::flex_double fraction(array_size);

    // The range of data frames and grid points to iterate over
    int j0 = roi_z[0];
    int j1 = roi_z[1];
    int v30 = - grid_origin_e3_;
    int v31 = + grid_origin_e3_;
    
    // A constant used in the solution to the integrals below. 
    double sigr2 = 1.0 / (std::sqrt(2.0) * (sigma_mosaicity_ / 
                                            std::abs(zeta)));

    // Loop over all j data frames in the region around the reflection
    for (int i = 0, j = j0; j <= j1; ++j) {

        // The data frame j covers the range of phi such that
        // rj = {phi':phi0 + (j-1)dphi <= phi' >= phi0 + jdpi}
        // Therefore the range of phi for j is given as follows.
        double bj = starting_angle_ + (j + starting_frame_) * oscillation_range_;
        double aj = bj - oscillation_range_;
        
        // Calculate the integral over rj (leaving out scaling factors):
        // I[exp(-(phi' - phi)^2 / (2 sigma^2)]
        double integral_j = (boost::math::erf((bj - phi) * sigr2) - 
                             boost::math::erf((aj - phi) * sigr2));
        
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
                double bv3 = ((v3 - 0.5) * step_size_e3_) / zeta + phi;
                double av3 = ((v3 + 0.5) * step_size_e3_) / zeta + phi;
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
                    fraction[i] = (boost::math::erf((bv3j - phi) * sigr2) - 
                                   boost::math::erf((av3j - phi) * sigr2)) * 
								   integral_j_r;
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
class XdsTransformDetectorBeamVectors {

public:

    /** The default constructor */
    XdsTransformDetectorBeamVectors() {}    
    
    /**
     * Initialise the class.
     * @param detector The detector struct
     * @param wavelength The wavelength of the radiation
     * @param n_div The number of pixel sub-divisions
     */
    XdsTransformDetectorBeamVectors(const equipment::Detector &detector,
                                    double wavelength,
                                    int n_div)
        : size_(detector.get_size()),
          origin_(detector.get_origin()),
          x_axis_(detector.get_x_axis().normalize() * 
                detector.get_pixel_size()[0]),
          y_axis_(detector.get_y_axis().normalize() * 
                detector.get_pixel_size()[1]),
          normal_(detector.get_normal()),
          distance_(detector.get_distance()),
          wavelength_r_(1.0 / wavelength),
          n_div_(n_div) {}

    /**
     * Calculate the beam vector at every pixel on the detector, sub-divided
     * into (n_div * n_div) equal areas. This is done to remove a certain
     * amount of processing from being done per reflection and ensuring it
     * is only done before the reflections are procesed.
     * @returns An array of beam vectors 
     */
    af::flex_vec3_double calculate() {
    
        // Calculate the image size
        int x_size = size_[0] * n_div_;
        int y_size = size_[1] * n_div_;
        int image_size = x_size * y_size;
        double n_div_r = 1.0 / (double)n_div_;
        
        // Create the necessary arrays
        af::flex_vec3_double detector_s1x(x_size);
        af::flex_vec3_double detector_s1y(y_size);
        af::flex_vec3_double detector_s1(image_size);
        
        // Calculate the x and y components of the detector s1 vectors
        for (std::size_t i = 0; i < x_size; ++i) {
            detector_s1x[i] = (i * n_div_r - origin_[0]) * x_axis_;
        }
        for (std::size_t j = 0; j < y_size; ++j) {
            detector_s1y[j] = (j * n_div_r - origin_[1]) * y_axis_;
        }

        // Calculate the s1 vector for each sub-division of the detector        
        for (std::size_t k = 0, j = 0; j < y_size; ++j) {
            for (std::size_t i = 0; i < x_size; ++i, ++k) {
                detector_s1[k] = (
                    detector_s1x[i] + 
                    detector_s1y[j] + 
                    distance_ * normal_).normalize() * wavelength_r_;
            }
        }
        
        // Return the s1 vector
        return detector_s1;
    }

private:
    
    equipment::Detector detector;
    scitbx::vec2 <int> size_;
    scitbx::vec2 <double> origin_;
    scitbx::vec3 <double> x_axis_;
    scitbx::vec3 <double> y_axis_;
    scitbx::vec3 <double> normal_;
    double distance_;
    double wavelength_r_;
    int n_div_;
};

/**
 * Class representing the XDS transform of the reflection profile on the
 * detector to the XDS reciprocal lattice coordinate frame.
 */
class XdsTransform {

public:

    /** The default constructor */
    XdsTransform() {}

    /**
     * Initialise the transform.
     * @param grid The transform grid container
     * @param image The raw image volume
     * @param mask THe reflection mask
     * @param detector The detector struct
     * @param beam The beam struct
     * @param gonio The goniometer struct
     * @param n_div The number of pixel sub divisions to use (default 5)
     */
    XdsTransform(XdsTransformGrid &grid,
                 const scitbx::af::flex_int &image,
                 const scitbx::af::flex_int &mask,
                 const equipment::Detector &detector,
                 const equipment::Beam &beam,
                 const equipment::Goniometer &gonio,
                 int n_div = 5)
        : e3_fraction_(gonio, grid),
          detector_s1_(
            XdsTransformDetectorBeamVectors(
                detector, 
                beam.get_wavelength(), 
                n_div).calculate()),
          grid_(grid),
          image_(image),
          mask_(mask),
          grid_size_(grid.get_size()),
          grid_origin_(grid.get_origin()),
          step_size_(grid.get_step_size()),
          s0_(beam.get_direction()),
          m2_(gonio.get_rotation_axis()),
          n_div_(n_div)
    {
        // Check the input
        DIALS_ASSERT(n_div > 0);
        DIALS_ASSERT(image.accessor().all().size() == 3);
        
        // Set image size
        image_size_ = scitbx::vec3 <int> (
            image.accessor().all()[0],
            image.accessor().all()[1],
            image.accessor().all()[2]);
    }

    void calculate(int reflection_index,
                   scitbx::af::tiny <int, 6> roi, 
                   scitbx::vec3 <double> s1,
                   double phi);
    
    /**
     * Calculate the XDS transform for all reflections
     * @param xyz The image volume coordinates
     * @param s1 The diffracted beam vectors
     * @param phi The rotation angles
     */  
    void calculate(const af::flex_tiny6_int &roi, 
                   const af::flex_vec3_double &s1, 
                   const scitbx::af::flex_double &phi) {
        DIALS_ASSERT(roi.size() == s1.size());
        DIALS_ASSERT(roi.size() == phi.size());
        DIALS_ASSERT(roi.size() <= grid_.get_n_reflections());
        for (int i = 0; i < roi.size(); ++i) {
            calculate(i, roi[i], s1[i], phi[i]);
        }
    }

private:

    XdsTransformE3Fraction e3_fraction_;
    af::flex_vec3_double detector_s1_;
    XdsTransformGrid grid_;
    scitbx::af::flex_int image_;
    scitbx::af::flex_int mask_;
    scitbx::vec3 <int> image_size_;
    scitbx::vec3 <int> roi_size_;
    scitbx::vec3 <int> grid_size_;
    scitbx::vec3 <int> grid_origin_;
    scitbx::vec3 <double> step_size_;
    scitbx::vec3 <double> s0_;
    scitbx::vec3 <double> m2_;
    int n_div_;
};

/**
 * Transform the profile of the reflection at detector point xyz, with beam
 * vector s1 and rotation angle phi, to the XDS reciprocal lattice coordinate
 * frame.
 *
 * We treat the image pixels as bins of a histogram. In effect we want to 
 * redistribute the image pixel counts to the XDS grid elements. We do this
 * by transforming the detector pixel coordinates around the reflection to the
 * XDS reciprocal lattice coordinate frame and determining the fraction of the
 * pixel value that is given to each element in the XDS grid.
 *
 * @param reflection_index The index of the reflection
 * @param xyz The coordinate of the reflection in the detector image volume
 * @param s1 The beam vector of the reflection
 * @param phi The rotation angle of the reflection
 * @throws std::rumtime_error if input is invalid
 */
void XdsTransform::calculate(int reflection_index,
                             scitbx::af::tiny <int, 6> roi,
                             scitbx::vec3 <double> s1,
                             double phi)
{
    // Check the reflection index
    DIALS_ASSERT(reflection_index >= 0 && 
                 reflection_index < grid_.get_n_reflections());
    
    // Constant for scaling values
    static const double r2d = 1.0 / scitbx::constants::pi_180;

    // Calculate the strides for indexing multidimensional arrays
    int div_image_stride_x = image_size_[2] * n_div_;
    int image_stride_x = image_size_[2];
    int image_stride_y = image_size_[1] * image_stride_x;
    int grid_stride_c1 = grid_size_[2];
    int grid_stride_c2 = grid_size_[1] * grid_stride_c1;
    int grid_stride_c3 = grid_size_[0] * grid_stride_c2;
    int grid_offset = reflection_index * grid_stride_c3;

    // Get the grid data array
    scitbx::af::flex_double grid = grid_.get_data();

    // Calculate the x, y, z ranges to iterate over
    int x0 = roi[0] * n_div_;
    int x1 = roi[1] * n_div_;
    int y0 = roi[2] * n_div_;
    int y1 = roi[3] * n_div_;
    int z0 = roi[4];
    int z1 = roi[5];
    
    // Check the data range
    DIALS_ASSERT(x0 >= 0 && x1 < image_size_[2] * n_div_ &&
                 y0 >= 0 && y1 < image_size_[1] * n_div_ &&
                 z0 >= 0 && z1 < image_size_[0]);
    
    // Calculate 1 / n_div and 1 / (n_div*n_div) for convenience
    double n_div_r = 1.0 / n_div_;
    double div_fraction = n_div_r * n_div_r;

    // Calculate the reflection coordinate system e1 and e2 axes, and zeta, the
    // loretz correction (used to calculate component on e3 axis
    scitbx::vec3 <double> e1 = s1.cross(s0_).normalize();
    scitbx::vec3 <double> e2 = s1.cross(e1).normalize();
    double zeta = m2_ * e1;
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
    scitbx::af::flex_double fraction = e3_fraction_.calculate( 
        scitbx::vec2 <int> (z0, z1), phi, zeta); 
    
    // Loop through all the pixels (and their sub-divisions). Calculate the
    // coordinate of each pixel in the XDS coordinate frame e1 and e2 axes.
    // Find the grid point in which the calculate point is contained and then
    // add the counts for that pixel to the grid. See Kabsch 2010
    for (int yy = y0; yy < y1; ++yy) {
        for (int xx = x0; xx < x1; ++xx) {
            double c1 = e1 * detector_s1_[xx + yy * div_image_stride_x] - c11;
            double c2 = e2 * detector_s1_[xx + yy * div_image_stride_x] - c21;
            int gi = grid_origin_[2] + c1 / step_size_[2];
            int gj = grid_origin_[1] + c2 / step_size_[1];
            if (gi < 0 || gi >= 9 || gj < 0 || gj >= 9) {
                continue;
            }
            int x = xx * n_div_r;
            int y = yy * n_div_r;
            int image_index = x + y * image_stride_x + z0 * image_stride_y;
            int fraction_index = 0;
            for (int z = z0; z <= z1; ++z) {
                if (mask_[image_index] == reflection_index) {
                    int value = image_[image_index] * div_fraction;
                    int grid_index = grid_offset + gi + gj * grid_stride_c1;
                    for (int gk = 0; gk < grid_size_[0]; ++gk) {
                        grid[grid_index] += value * fraction[fraction_index];
                        grid_index += grid_stride_c2;
                        fraction_index++;
                    }
                } else {
                    fraction_index += grid_size_[0];
                }
                image_index += image_stride_y;
            }
        }
    }
}

}} // namespace = dials::integration

#endif // DIALS_INTEGRATION_XDS_TRANSFORM_H
