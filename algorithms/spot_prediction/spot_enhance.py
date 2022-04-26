from __future__ import annotations

from collections import Counter

import numpy as np
import numpy.ma as ma
import plotly.express as px
from scipy.spatial.transform import Rotation as R

# pixel depth
D = 1  # mm

# pixel size
px_x = 0.172  # mm
px_y = 0.172  # mm

# detector normal
det_normal = np.array([-0.00186461, 0.00883524, -0.99995923])

# absorption factor
mu = 0.6245476712834831  # mm^-1
# absorption length
s = 1 / mu


def lab_to_detector(vector):
    """
    From lab to detector there is a pi rotation around x
    plus a slight detector tilt
    (see Jur van den Berg (https://math.stackexchange.com/users/91768/jur-van-den-berg),
    Calculate Rotation Matrix to align Vector A to Vector B in 3d?,
    URL (version: 2016-09-01): https://math.stackexchange.com/q/476311)
    """
    r_flipped = R.from_rotvec(np.pi * np.array([1, 0, 0]))

    z_lab = np.array([0, 0, 1])
    v = np.cross(z_lab, det_normal)
    v_mat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])

    c = z_lab.dot(det_normal)

    r = R.from_matrix(np.eye(3) + v_mat + v_mat * v_mat / (1 + c))

    return r_flipped.apply(r.apply(vector))


def map_mm_to_px(position_array, px_x, px_y):
    """
    position_array : ndarray, size=(3,n)
    """
    index_x = ((position_array.T[0]) // px_x).astype(int)
    index_y = ((position_array.T[1]) // px_y).astype(int)

    pixels = list(zip(index_x, index_y))
    return pixels


def match_to_px_grid(array, s1_map):
    return np.repeat(np.repeat(array, s1_map, axis=0), s1_map, axis=1)


def position_grid(num_pix, px_x, px_y, s1_map):
    """
    make a grid for pixel positions which will be used as s1 origin.
    If s1_map is 1 then the pixel position is its center
    Else, grid s1_map x s1_map every pixel
    """
    px_pos_x, px_pos_y = np.meshgrid(
        np.arange(px_x / (2 * s1_map), px_x * num_pix, px_x / s1_map),
        np.arange(px_y / (2 * s1_map), px_y * num_pix, px_y / s1_map),
    )

    px_pos_z = np.zeros((num_pix * s1_map, num_pix * s1_map))

    # numpy arrays have axis y before x: v[y,x]
    px_pos = np.stack((px_pos_y, px_pos_x, px_pos_z)).T
    return px_pos


def plot_view_in_det_frame(intensity_2D_array, pad, title):
    fig = px.imshow(intensity_2D_array, title=title, zmin=0)
    # add padding box
    pad_edge = pad - 0.5
    fig.add_scatter(
        x=[pad_edge, pad_edge + 7, pad_edge + 7, pad_edge, pad_edge],
        y=[pad_edge, pad_edge, pad_edge + 7, pad_edge + 7, pad_edge],
        mode="lines",
    )
    fig.show()


def pad_2Darray(array, pad):
    return np.pad(array, ((pad, pad), (pad, pad)), "constant", constant_values=0)


def pad_3Darray(array, pad):
    """
    pad every slice in 3D array
    """
    return np.pad(
        array, ((pad, pad), (pad, pad), (0, 0)), "constant", constant_values=0
    )


def psf(s1_vectors, intensity, px_size, padding, num_MC, s1_map=1):
    int_padded = pad_2Darray(intensity, padding)
    s1_padded = pad_3Darray(s1_vectors, padding)

    num_pix = len(int_padded)
    px_pos = position_grid(num_pix, *px_size, s1_map)

    # if s1_map>1 repeat the s1 and int array such that they match px_pos
    int_matched = match_to_px_grid(int_padded, s1_map)
    s1_matched = match_to_px_grid(s1_padded, s1_map)

    # ravel all array to shape = (padded_area_length + padded_area_height, quantity_dimensions)
    int_flat = int_matched.ravel()
    total_length = len(int_flat)
    s1_flat = s1_matched.reshape(total_length, 3)
    px_pos_flat = px_pos.reshape(total_length, 3)

    # convert s1 vectors from lab to detector frame
    s1_det_flat = lab_to_detector(s1_flat, det_normal)

    # mask zeros in spot region
    int_nonzero = ma.masked_values(int_flat, 0)

    summed_int = int(np.sum(int_nonzero))
    # summed_int x num_MC array of random numbers
    rnd = np.random.rand(summed_int * num_MC)

    # summed_int x num_MC array of random paths traveled
    l = -s * np.log(rnd)

    # keep data just in non-zero region
    int_valid_region = int_nonzero[~int_nonzero.mask]
    s1_valid_region = s1_det_flat[~int_nonzero.mask]
    pos_valid_region = px_pos_flat[~int_nonzero.mask]

    s1_per_hit = np.repeat(s1_valid_region, int_valid_region, axis=0)
    s1_per_MC = np.repeat(s1_per_hit, num_MC, axis=0)

    final_pos = s1_per_MC.T * l

    pos_per_hit = np.repeat(pos_valid_region, int_valid_region, axis=0)
    pos_per_MC = np.repeat(pos_per_hit, num_MC, axis=0)
    coords = pos_per_MC + final_pos.T

    # ignore x-rays that escape through the back of detector
    escaping_zplus = ma.masked_greater(coords[:, 2], D)

    # ignore x-rays that travel outside padded region
    escaping_xplus = ma.masked_greater(coords[:, 0], px_size[0] * num_pix)
    escaping_xminus = ma.masked_less_equal(coords[:, 0], 0)
    escaping_yplus = ma.masked_greater(coords[:, 1], px_size[1] * num_pix)
    escaping_yminus = ma.masked_less_equal(coords[:, 1], 0)

    # if escaped mask away
    escaping_mask = (
        escaping_xminus.mask
        | escaping_xplus.mask
        | escaping_yminus.mask
        | escaping_yplus.mask
        | escaping_zplus.mask
    )

    indices = map_mm_to_px(coords[~escaping_mask], *px_size)

    # count hits per pixel
    c = Counter(indices)
    key_x, key_y = zip(*c.keys())

    counts_on_back = np.zeros((num_pix, num_pix))
    # populate count array with counts
    counts_on_back[[key_x], [key_y]] = list(c.values())

    # counts_on_back[counts_on_back==0] =np.nan

    # plot some things
    plot_view_in_det_frame(
        int_padded.reshape(num_pix, num_pix), padding, "on the front of the detector"
    )
    plot_view_in_det_frame(
        counts_on_back.T / (num_MC * s1_map * s1_map),
        padding,
        "at the back of the detector",
    )

    # plot_vector_in_det_frame(s1_valid_region)
    # plot_path_length_dist(final_pos[2][final_pos[2]<=D])
    det_efficiency = np.sum(counts_on_back) / summed_int / (num_MC)

    return counts_on_back / (num_MC * s1_map * s1_map), det_efficiency


def enhance(spot, s1, pixel_param):
    """
    Given a spot as measured on the back of the detector
    deconvolve the detector effect to recover the spot
    on the front of the detector.

    :param spot: array
            2D array of intensity values contained in a shoebox
            as measured on the 'back' of the detector
    :param s1: array
            3D array: 2D array of incident array vectors
    :param pixel_param: tuple
            pixel parameters including size and thickness, DQE and absorption length

    :return enhanced_spot: array
    """

    # for the deconvolution to recover exactly the spot
    # on the front of the detector the back of the detector
