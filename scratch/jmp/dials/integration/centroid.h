
#ifndef DIALS_INTEGRATION_CENTROID_H
#define DIALS_INTEGRATION_CENTROID_H

#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include <scitbx/array_family/tiny.h>
#include "../error.h"

namespace dials { namespace integration {

/** Calculate the centroid of a 2D image */
template <typename ImageType>
scitbx::vec2 <double> 
centroid2d(const ImageType &image, scitbx::af::tiny <int, 4> roi) {
    double xc = 0, yc = 0, count = 0;
    for (std::size_t j = roi[2]; j < roi[3]; ++j) {
        for (std::size_t i = roi[0]; i < roi[1]; ++i) {
            xc += (i+0.5) * image(j, i);
            yc += (j+0.5) * image(j, i);
            count += image(j, i);
        }
    }
    DIALS_ASSERT(count > 0);    
    return scitbx::vec2 <double> (xc / count, yc / count);
}

/** Calculate the centroid of a 3D image */
template <typename VolumeType>
scitbx::vec3 <double> 
centroid3d(const VolumeType &image, scitbx::af::tiny <int, 6> roi) {
    double xc = 0, yc = 0, zc = 0, count = 0;
    for (std::size_t k = roi[4]; k < roi[5]; ++k) {
        for (std::size_t j = roi[2]; j < roi[3]; ++j) {
            for (std::size_t i = roi[0]; i < roi[1]; ++i) {
                xc += (i+0.5) * image(k, j, i);
                yc += (j+0.5) * image(k, j, i);
                zc += (k+0.5) * image(k, j, i);
                count += image(k, j, i);
            }
        }
    }
    DIALS_ASSERT(count > 0);       
    return scitbx::vec3 <double> (xc / count, yc / count, zc / count);
}

/** Calculate the centroid of a 2D image with a mask */
template <typename ImageType, typename MaskType>
scitbx::vec2 <double> 
centroid2d(const ImageType &image, const MaskType &mask, 
         scitbx::af::tiny <int, 4> roi, int value) {
    double xc = 0, yc = 0, count = 0;
    for (std::size_t j = roi[2]; j < roi[3]; ++j) {
        for (std::size_t i = roi[0]; i < roi[1]; ++i) {
            if (mask(j, i) == value) {
                xc += (i+0.5) * image(j, i);
                yc += (j+0.5) * image(j, i);
                count += image(j, i);
            }
        }
    }
    DIALS_ASSERT(count > 0);    
    return scitbx::vec2 <double> (xc / count, yc / count);
}

/** Calculate the centroid of a 3D image with a mask */
template <typename VolumeType, typename MaskType>
scitbx::vec3 <double> 
centroid3d(const VolumeType &image, const MaskType &mask, 
         scitbx::af::tiny <int, 6> roi, int value) {
    double xc = 0, yc = 0, zc = 0, count = 0;
    for (std::size_t k = roi[4]; k < roi[5]; ++k) {
        for (std::size_t j = roi[2]; j < roi[3]; ++j) {
            for (std::size_t i = roi[0]; i < roi[1]; ++i) {
                if (mask(k, j, i) == value) {
                    xc += (i+0.5) * image(k, j, i);
                    yc += (j+0.5) * image(k, j, i);
                    zc += (k+0.5) * image(k, j, i);
                    count += image(k, j, i);
                }
            }
        }
    }
    DIALS_ASSERT(count > 0);
    return scitbx::vec3 <double> (xc / count, yc / count, zc / count);
}

}} // namespace dials::integration

#endif // DIALS_INTEGRATION_CENTROID_H
