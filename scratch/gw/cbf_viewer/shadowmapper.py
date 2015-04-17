# -*- coding: utf-8 -*-

""" Contains the Shadow class. This class is used to calculate the
mapping of the complex goniometer shadow onto the detector. The shadow
is drawn on the diffraction image texture using PIL.
"""

import numpy as np
from numpy.linalg import norm

import PIL.ImageDraw
import PIL.ImageChops

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import dxtbx

import analyser
import viewer

class Shadow(object):

    def  __init__(self, gonio_coords, filename, detector_angle,
                  detector_distance, origin, fast, slow,
                  offset_matrix = np.zeros(3)):

        self.gonio_coords = gonio_coords
        self.filename = filename
        self.detector_angle = detector_angle
        self.detector_distance = detector_distance
        self.origin = origin
        self.fast = fast
        self.slow = slow
        self.offset_matrix = offset_matrix

    def orient_gonio(self, kappa = -180.0, omega = 90.0):
        # This code makes the angles the right way round - this is
        # necessary because the texture is flipped later on by OpenGL.
        newkappa = -180.0 + (-180.0 - kappa)
        newomega = 90.0 + (90.0 - omega)

        self.k = np.deg2rad(newkappa)
        self.om = np.deg2rad(newomega)
        u = np.cos(np.deg2rad(50.0))
        v = np.sin(np.deg2rad(50.0))

        sink = np.sin(self.k)
        cosk = np.cos(self.k)
        cosom = np.cos(self.om)
        sinom = np.sin(self.om)
        # r_k is the kappa rotation matrix.
        r_k = np.array([[u ** 2 + ((v ** 2) * cosk), u * v * (1.0 - cosk),
                        -v * sink], [u * v * (1 - cosk),
                        v ** 2 + ((u ** 2) * cosk), u * sink],
                        [v * sink, -u * sink, (u ** 2 + v ** 2) * cosk]])
        # rot_k is the gonio coords after kappa rotation.
        rot_k = np.dot(self.gonio_coords, r_k)

        for i in range(len(rot_k)):
            for j in range(3):
                rot_k[i, j] += self.offset_matrix[j]
        # r_om is the omega rotation matrix.
        r_om = np.array([[1.0, 0.0, 0.0],
                        [0.0, cosom, sinom],
                        [0.0, -sinom, cosom]])

        new_coords = np.dot(rot_k, r_om)
        return new_coords

    def shadow_texture(self, kappa = -180.0, omega = 90.0, blank = False):
        '''Extracts the diffraction image from the cbf file, then draws
        the shadow projection on it using PIL.
        The blank arg is used to toggle a blank detector with shadow.
        '''

        myimage = analyser.Image(self.filename)
        pix, res = myimage.detector_parameters()[4 : 7 : 2]
        # The code below converts the diffraction image to a numpy
        # array unless blank is true, in which case a blank detector
        # is displayed with shadow.
        if blank is False:
            img_array = myimage.image_to_array()

        else:
            # The blank black image initialised below is called white
            # because it will later be inverted.
            white = PIL.Image.new("L", (res), "black")
            img_array = np.asarray(white)

        tex = PIL.Image.fromarray(img_array)
        tex_shadow = PIL.Image.fromarray(img_array)

        shadow_coords = ()
        normal = np.cross(self.fast[0], self.slow[0])
        normal /= norm(normal)

        new_coords = self.orient_gonio(kappa, omega)
        normal = np.cross(self.fast[0], self.slow[0])
        normal /= norm(normal)
        dist = np.zeros(len(new_coords))
        proj = np.ndarray(np.shape(new_coords))
        x_det2 = np.zeros(len(new_coords))
        y_det2 = np.zeros(len(new_coords))

        for i in range(len(dist)):
            dist[i] = (np.dot(self.origin[0], normal) /
                    np.dot(new_coords[i], normal))
            proj[i] = dist[i] * new_coords[i]

            x_det2[i] = (-np.sign(self.origin[0, 0] - proj[i, 0]) *
                       (norm(np.cross(self.slow[0], (self.origin[0]
                        - proj[i]))) / norm(self.slow[0])) / pix[0])
            y_det2[i] = (np.sign(self.origin[0, 1] - proj[i, 1]) *
                       (norm(np.cross(self.fast[0], (self.origin[0]
                        - proj[i]))) / norm(self.fast[0])) / pix[1])

            # This line prevents the entire detector being covered by
            # the shadow when it shouldn't be. This occurs when a z
            # coord in new_coords is very close to zero.
            if dist[i] < 0.0 and abs(new_coords[i, 2]) > 10.0 ** -10:
                shadow_coords += (x_det2[i], y_det2[i])

        draw = PIL.ImageDraw.Draw(tex_shadow)
        if len(shadow_coords) > 2:
            draw.polygon(shadow_coords, fill = 255)

        tex1 = tex_shadow.resize((512, 512))
        tex2 = PIL.ImageChops.invert(tex1)

        tex_data = np.array(list(tex2.getdata()), np.uint8)
        tex_noshad = PIL.ImageChops.invert(tex.resize((512, 512)))
        noshad_data = np.array(list(tex_noshad.getdata()), np.uint8)

        return tex_data, tex, tex2, tex_shadow, tex_noshad, noshad_data

    def affected_pixels(self, kappa = -180.0, omega = 90.0, blank = False):
        '''Counts pixels affected by goniometer shadow. Returns the
        difference between np arrays of shadowed and unshadowed
        textures, the number of affected pixels, and the position of
        each affected pixel.
        '''
        (tex_data, tex, tex2, tex_shadow, tex_noshad,
        noshad_data) = self.shadow_texture(kappa, omega, blank)

        tex_array = np.asarray(tex)
        tex_shadow_array = np.asarray(tex_shadow)

        diffs = tex_shadow_array - tex_array
        pixels = np.count_nonzero(diffs)
        indices = np.nonzero(diffs)

        return diffs, pixels, indices

    def pixel_plotter(self, blank = True, num_points = 50):
        ''' This function plots a 3D graph of percentage of affected
        pixels as a function of omega and kappa. This takes a long
        time.
        '''
        kapparange = np.linspace(-190.0, 10.0, num_points)
        omegarange = np.linspace(-270.0, 270.0, num_points)
        X, Y = np.meshgrid(kapparange, omegarange)
        pixels = (np.array([self.affected_pixels(x, y, blank=blank)[1]
        for x, y in zip(np.ravel(X), np.ravel(Y))]) /
                    (viewer.dimensions[0][0] * viewer.dimensions[0][1] /
                    (viewer.pix[0] * viewer.pix[1] * 100.0)))

        Z = pixels.reshape(X.shape)

        ax = plt.contourf(X, Y, Z)
        cbar = plt.colorbar()
        cbar.set_label('Percentage of detector in shadow')

        plt.show()
        plt.savefig('pixplot.png', dpi=96, transparent=True)

day2_coords = np.array([[107.61, 0.00, -40.00],
[107.61, -6.95, -39.39],
[107.61, -13.68, -37.59],
[107.61, -20.00, -34.64],
[107.61, -25.71, -30.64],
[107.61, -30.64, -25.71],
[107.61, -34.64, -20.00],
[107.61, -37.59, -13.68],
[107.61, -39.39, -6.95],
[107.61, -40.00, 0.00],
[107.61, -39.39, 6.95],
[107.61, -37.59, 13.68],
[107.61, -34.64, 20.00],
[107.61, -30.64, 25.71],
[107.61, -25.71, 30.64],
[107.61, -20.00, 34.64],
[107.61, -13.68, 37.59],
[107.61, -6.95, 39.39],
[107.61, 0.00, 40.00],
[107.49, 7.84, 39.49],
[107.39, 15.69, 38.97],
[107.27, 23.53, 38.46],
[107.16, 31.37, 37.94],
[101.76, 33.99, 36.25],
[96.37, 36.63, 34.56],
[90.98, 39.25, 33.00],
[85.58, 41.88, 31.18],
[80.89, 47.06, 31.00],
[76.55, 51.51, 31.03],
[72.90, 55.04, 31.18],
[66.86, 60.46, 31.67],
[62.10, 64.41, 32.25],
[57.22, 68.19, 33.00],
[52.83, 71.88, 32.50],
[48.57, 75.45, 31.01],
[44.58, 78.80, 28.58],
[40.97, 81.83, 25.28],
[37.86, 84.44, 21.21],
[35.33, 86.56, 16.50],
[33.47, 88.13, 11.29],
[32.33, 89.08, 5.73],
[31.94, 89.41, 0.00],
[32.33, 89.08, -5.73],
[33.47, 88.13, -11.29],
[35.33, 86.56, -16.50],
[37.86, 84.44, -21.21],
[40.97, 81.83, -25.28],
[44.58, 78.80, -28.58],
[48.57, 75.45, -31.01],
[52.83, 71.88, -32.50],
[57.22, 68.19, -33.00],
[62.10, 64.41, -32.25],
[66.86, 60.46, -31.67],
[72.90, 55.04, -31.18],
[76.55, 51.51, -31.03],
[80.89, 47.06, -31.00],
[85.58, 41.88, -31.18],
[90.98, 39.25, -33.00],
[96.37, 36.63, -34.56],
[101.76, 33.99, -36.25],
[107.16, 31.37, -37.94],
[107.27, 23.53, -38.46],
[107.39, 15.69, -38.97],
[107.49, 7.84, -39.49],
[107.61, 0.00, -40.00]])

if __name__ == "__main__":
    myshadow = Shadow(day2_coords, 250.0, -180.0, 90.0)
    print myshadow.get_projection()
