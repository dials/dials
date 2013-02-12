
#ifndef DIALS_INTEGRATION_CENTROID_H
#define DIALS_INTEGRATION_CENTROID_H

#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include <scitbx/array_family/tiny.h>
#include <scitbx/array_family/flex_types.h>
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

template <typename MaskType>
void mask_adjacent_pixels(MaskType &mask, scitbx::vec3 <std::size_t> index) {
    typedef scitbx::vec3 <std::size_t> index_type;
    DIALS_ASSERT(mask.accessor().all().size() == 3);
    scitbx::af::small <long int, 10> size = mask.accessor().all();
    int i = index[0], j = index[1], k = index[2];
    if (i >= 0 && i < size[2] && j >= 0 && j < size[1] && k >= 0 && k < size[0]) {
        if (mask(k, j, i) == 1) {
            mask(k, j, i) = 2;
            mask_adjacent_pixels(mask, index_type(i-1, j, k));
            mask_adjacent_pixels(mask, index_type(i+1, j, k));
            mask_adjacent_pixels(mask, index_type(i, j-1, k));
            mask_adjacent_pixels(mask, index_type(i, j+1, k));
            mask_adjacent_pixels(mask, index_type(i, j, k-1));
            mask_adjacent_pixels(mask, index_type(i, j, k+1));
        }
    }
}

template <typename VolumeType, typename MaskType>
scitbx::vec3 <double>
centroid_reflection(const VolumeType &image, const MaskType &mask,
                    scitbx::af::tiny <int, 6> roi, int value, double threshold) {
    typedef scitbx::vec3 <std::size_t> index_type;

    scitbx::af::flex_int centroid_mask(scitbx::af::flex_grid <> (
                                            roi[5]-roi[4],
                                            roi[3]-roi[2],
                                            roi[1]-roi[0]), 0);

    scitbx::af::flex_int spot_image(scitbx::af::flex_grid <> (
                                            roi[5]-roi[4],
                                            roi[3]-roi[2],
                                            roi[1]-roi[0]), 0);

    // Threshold above given level and find the maximum value
    index_type max_index(-1, -1, -1);
    double max_image = 0;
    for (std::size_t k = roi[4]; k < roi[5]; ++k) {
        for (std::size_t j = roi[2]; j < roi[3]; ++j) {
            for (std::size_t i = roi[0]; i < roi[1]; ++i) {
                if (mask(k, j, i) == value) {
                    if (image(k, j, i) >= threshold) {
                        std::size_t k1 = k - roi[4];
                        std::size_t j1 = j - roi[2];
                        std::size_t i1 = i - roi[0];
                        centroid_mask(k1, j1, i1) = 1;
                        spot_image(k1, j1, i1) = image(k, j, i);
                        if (max_image < image(k, j, i)) {
                            max_image = image(k, j, i);
                            max_index = index_type(i1, j1, k1);
                        }
                    }
                }
            }
        }
    }

    // Ensure we have an index
    DIALS_ASSERT(max_image > 0);

    // Select the spot peak pixels
    mask_adjacent_pixels(centroid_mask, max_index);

    // Create the roi
    scitbx::af::tiny <int, 6> spot_roi(0, roi[1]-roi[0], 0, roi[3]-roi[2], 0, roi[5]-roi[4]);

    // Centroid the spot image
    return centroid3d(spot_image, centroid_mask, spot_roi, 2) +
        scitbx::vec3 <double> (roi[0], roi[2], roi[4]);
}

template <typename VolumeType>
scitbx::vec3 <double>
centroid_reflection(const VolumeType &image, scitbx::af::tiny <int, 6> roi,
                    double threshold) {
    typedef scitbx::vec3 <std::size_t> index_type;

    scitbx::af::flex_int centroid_mask(scitbx::af::flex_grid <> (
                                            roi[5]-roi[4],
                                            roi[3]-roi[2],
                                            roi[1]-roi[0]), 0);

    scitbx::af::flex_int spot_image(scitbx::af::flex_grid <> (
                                            roi[5]-roi[4],
                                            roi[3]-roi[2],
                                            roi[1]-roi[0]), 0);

    // Threshold above given level and find the maximum value
    index_type max_index(-1, -1, -1);
    double max_image = 0;
    for (std::size_t k = roi[4]; k < roi[5]; ++k) {
        for (std::size_t j = roi[2]; j < roi[3]; ++j) {
            for (std::size_t i = roi[0]; i < roi[1]; ++i) {
                if (image(k, j, i) >= threshold) {
                    std::size_t k1 = k - roi[4];
                    std::size_t j1 = j - roi[2];
                    std::size_t i1 = i - roi[0];
                    centroid_mask(k1, j1, i1) = 1;
                    spot_image(k1, j1, i1) = image(k, j, i);
                    if (max_image < image(k, j, i)) {
                        max_image = image(k, j, i);
                        max_index = index_type(i1, j1, k1);
                    }
                }
            }
        }
    }

    // Ensure we have an index
    DIALS_ASSERT(max_image > 0);

    // Select the spot peak pixels
    mask_adjacent_pixels(centroid_mask, max_index);

    // Create the roi
    scitbx::af::tiny <int, 6> spot_roi(0, roi[1]-roi[0], 0, roi[3]-roi[2], 0, roi[5]-roi[4]);

    // Centroid the spot image
    return centroid3d(spot_image, centroid_mask, spot_roi, 2) +
        scitbx::vec3 <double> (roi[0], roi[2], roi[4]);
}


}} // namespace dials::integration

#endif // DIALS_INTEGRATION_CENTROID_H
