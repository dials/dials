# -*- coding: utf-8 -*-

"""
This module is used to import and use the other four CBFviewer
modules. It contains only one function, which runs the program.
"""

import wx
import numpy as np

import analyser
import modeller
import viewer
import shadowmapper

det = []

def view_experiment(filename, sample_x_length=4.0, sample_y_length=0.5,
    sample_z_length=8.0):

    """
    Read a cbf image and display a 3D model of the experiment.

    This function is intended to be used by the block of code at the
    end of the module. This will start the program with a dialog box
    allowing the user to select a file.

    Arguments:
    filename -- the name of the cbf image to be read. Must include the
    .cbf extension.
    sample_x,y,z_length -- the dimensions of the sample in mm.

    The model can be rotated by clicking and dragging the mouse, and
    zoomed with the mousewheel. The model is displayed in a wxPython
    based OpenGLCanvas.

    Note that the default sample dimensions are much larger than the
    true dimensions; however for the purposes of the program it is best
    to have a crystal large enough to be visible.
    """

    myimage = analyser.Image(filename)

    viewer.B = myimage.find_beam_direction()

    (viewer.detector_origin, viewer.fast_axis, viewer.slow_axis,
    viewer.dimensions, viewer.pix, viewer.d,
    res) = myimage.detector_parameters()

    ((i, j, k), (alpha, beta, gamma),
    viewer.R) = myimage.get_sample_orientation()

    viewer.sample_angles = (alpha, beta, gamma)

    viewer.img_array = myimage.image_to_array()

    mysample = modeller.Sample(sample_x_length, sample_y_length,
                                  sample_z_length)

    viewer.r_f = mysample.rotate(i, j, k)

    for i in range(viewer.num_panels):
        det.append(modeller.Panel(viewer.dimensions[i, 0],
                                    viewer.dimensions[i, 1]))

        viewer.d_r_f[i] = (det[i].get_coords(viewer.detector_origin[i],
                            viewer.fast_axis[i], viewer.slow_axis[i]))

    myshadow = shadowmapper.Shadow(shadowmapper.day2_coords, filename,
    viewer.DET_ANGLE, viewer.d, viewer.detector_origin, viewer.fast_axis,
    viewer.slow_axis)

    viewer.filename, viewer.shadowsetting, viewer.blanksetting = (filename,
                                                                False, False)

    (viewer.kappa, viewer.omega, viewer.offset_matrix,
    viewer.sample_x) = -180.0, 90.0, np.zeros(3), 0.0

    (viewer.tex_data, viewer.tex, viewer.tex2, viewer.tex_shadow,
    viewer.tex_noshad, viewer.noshad_data) = myshadow.shadow_texture()

    viewer.app = wx.App()
    viewer.frame = viewer.MainWindow()
    viewer.app.MainLoop()

if __name__ == '__main__':

    app = wx.App()
    done = False
    wildcard = "cbf images (*.cbf)|*.cbf|All Files (*.*)|*.*"
    dlg = wx.FileDialog(message = "Select a cbf image", parent = None,
                wildcard = wildcard, defaultFile = "", style = wx.OPEN)
    # Making the user repeatedly try to pick a new file if their
    # selection is invalid, without closing the dialog until 'cancel'
    # is clicked.
    while not done:
        try:
            if dlg.ShowModal() == wx.ID_OK:
                filename = dlg.GetPath()
                view_experiment(filename)
                done = True
            elif dlg.ShowModal() == wx.ID_CANCEL:
                done = True
                dlg.Destroy()

        except IOError:
            message = "Invalid file. Please choose another."
            caption = "Error"
            wx.MessageBox(message, caption, wx.ICON_ERROR)

    app.MainLoop()
