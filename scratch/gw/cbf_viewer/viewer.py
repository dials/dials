# -*- coding: utf-8 -*-

"""
This module takes input coordinate arrays and draws a 3D image of
the experimental geometry. This is done using wxPython and an OpenGL
GLCanvas.
"""

import sys
from copy import deepcopy

import PIL.Image
import PIL.ImageDraw
import PIL.ImageEnhance
import PIL.ImageChops

import numpy as np
from numpy.linalg import norm

import wx
from wx import glcanvas
import wx.lib.agw.floatspin as fs

from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *

import modeller
import launcher
import shadowmapper

d_r_f = dict()

class MyCanvasBase(glcanvas.GLCanvas):

    def __init__(self, parent):
        glcanvas.GLCanvas.__init__(self, parent, -1)
        self.init = False
        # initial mouse position
        self.lastx = self.x = 30
        self.lasty = self.y = 30
        self.size = None
        self.Bind(wx.EVT_ERASE_BACKGROUND, self.OnEraseBackground)
        self.Bind(wx.EVT_SIZE, self.OnSize)
        self.Bind(wx.EVT_PAINT, self.OnPaint)
        self.Bind(wx.EVT_LEFT_DOWN, self.OnMouseDown)
        self.Bind(wx.EVT_LEFT_UP, self.OnMouseUp)
        self.Bind(wx.EVT_MOTION, self.OnMouseMotion)
        self.zf = 1.0           # zoomfactor

    def OnEraseBackground(self, event):
        pass # Do nothing, to avoid flashing on Windows.

    def OnSize(self, event):
        size = self.size = self.GetClientSize()
        self.aspect = size.width / size.height
        if self.GetContext():
            self.SetCurrent()
            glViewport(0, 0, size.width, size.height)
            # The following code prevents the model from stretching
            # when resized.
            glMatrixMode(GL_PROJECTION)
            glLoadIdentity()
            if size.width >= size.height:
                glOrtho(-100.0 * self.zf * self.aspect,
                         100.0 * self.zf * self.aspect,
                        -100.0 * self.zf, 100.0 * self.zf, -2000.0, 2000.0)
            else:
                glOrtho(-100.0 * self.zf, 100.0 * self.zf,
                        -100.0 * self.zf / self.aspect,
                         100.0 * self.zf / self.aspect, -2000.0, 2000.0)
        event.Skip()

    def OnPaint(self, event):
        dc = wx.PaintDC(self)
        self.SetCurrent()
        if not self.init:
            self.InitGL()
            self.init = True
        self.OnDraw()

    def OnMouseDown(self, evt):
        self.CaptureMouse()
        self.x, self.y = self.lastx, self.lasty = evt.GetPosition()

    def OnMouseUp(self, evt):
        self.ReleaseMouse()

    # Rotate view with mouse drag
    def OnMouseMotion(self, evt):
        if evt.Dragging() and evt.LeftIsDown():
            self.lastx, self.lasty = self.x, self.y
            self.x, self.y = evt.GetPosition()
            #print evt.GetPosition()
            #self.lastx = self.x
            #self.x = evt.GetPosition()[0]
            self.Refresh(False)

class ViewerCanvas(MyCanvasBase):
    def InitGL(self):
        # Set viewing projection
        self.Bind(wx.EVT_MOUSEWHEEL, self.OnZoom)
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        if self.aspect >=1.0:
            glOrtho(-100.0 * self.zf * self.aspect,
                     100.0 * self.zf * self.aspect,
                    -100.0 * self.zf, 100.0 * self.zf, -2000.0, 2000.0)
        else:
            glOrtho(-100.0 * self.zf, 100.0 * self.zf,
                    -100.0 * self.zf / self.aspect,
                     100.0 * self.zf / self.aspect, -2000.0, 2000.0)

        # position viewer
        glMatrixMode(GL_MODELVIEW)
        glTranslatef(-3.0, -3.0, -3.0)

        # position objects
        glRotatef(self.y, 1.0, 0.0, 0.0)
        glRotatef(self.x, 0.0, 1.0, 0.0)

        glLightfv(GL_LIGHT0, GL_AMBIENT, GLfloat_4(1.0, 1.0, 1.0, 1.0))
        glEnable(GL_DEPTH_TEST)
        glEnable(GL_LIGHTING)
        glEnable(GL_LIGHT0)

    def OnZoom(self, evt):
        if evt.GetWheelRotation() < 0:
            self.zf += 0.1
        # The 'and' statement below prevents the user from setting zf
        # to non-positive values, which was causing strange bugs.
        elif evt.GetWheelRotation() > 0 and self.zf > 0.2:
            self.zf -= 0.1

        # Now loading the projection with the new zoom factor
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        if self.aspect >=1.0:
            glOrtho(-100.0 * self.zf * self.aspect,
                     100.0 * self.zf * self.aspect,
                    -100.0 * self.zf, 100.0 * self.zf, -2000.0, 2000.0)
        else:
            glOrtho(-100.0 * self.zf, 100.0 * self.zf,
                    -100.0 * self.zf / self.aspect,
                     100.0 * self.zf / self.aspect, -2000.0, 2000.0)

        glMatrixMode(GL_MODELVIEW)
        glTranslatef(0.0, 0.0, 0.0)
        self.Refresh(False)

    def OnDraw(self):
        # clear color and depth buffers
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)

        # drawing the sample
        glBegin(GL_QUADS)
        glMaterial(GL_FRONT_AND_BACK, GL_AMBIENT, (0.2, 0.2, 0.2, 1.0))
        for i in range(6):          #for each face
            glNormal3f(r_f[i][0, 0], r_f[i][0, 1], r_f[i][0, 2])

            for j in range(1,5):    #for each vertex on the face
                glVertex3f(r_f[i][j, 0], r_f[i][j, 1], r_f[i][j, 2])
        glEnd()

        # Drawing the beam and rotation axis
        glBegin(GL_LINES)
        glVertex3f(0.0, 0.0, 0.0)
        glVertex3f(200.0 * B[0], 200.0 * B[1], 200.0 * B[2])

        glVertex3f(200.0 * R[0], 200.0 * R[1], 200.0 * R[2])
        glVertex3f(-200.0 * R[0], -200.0 * R[1], -200.0 * R[2])

        # Drawing the axes
        glVertex3f(0.0, 0.0, 0.0)
        glVertex3f(10.0, 0.0, 0.0)

        glVertex3f(0.0, 0.0, 0.0)
        glVertex3f(0.0, 10.0, 0.0)

        glVertex3f(0.0, 0.0, 0.0)
        glVertex3f(0.0, 0.0, 10.0)

        glEnd()

        # Drawing the axis arrows (this was a quick solution; probably
        # uses more lines of code than necessary).
        glBegin(GL_QUADS)
        glNormal3f(-1.0, 0.0, 0.0)
        glVertex3f(10.0, 1.0, 1.0)
        glVertex3f(10.0, 1.0, -1.0)
        glVertex3f(10.0, -1.0, -1.0)
        glVertex3f(10.0, -1.0, 1.0)

        glNormal3f(0.0, -1.0, 0.0)
        glVertex3f(1.0, 10.0, 1.0)
        glVertex3f(1.0, 10.0, -1.0)
        glVertex3f(-1.0, 10.0, -1.0)
        glVertex3f(-1.0, 10.0, 1.0)

        glNormal3f(0.0, 0.0, 10.0)
        glVertex3f(1.0, 1.0, 10.0)
        glVertex3f(-1.0, 1.0, 10.0)
        glVertex3f(-1.0, -1.0, 10.0)
        glVertex3f(1.0, -1.0, 10.0)

        glEnd()

        glBegin(GL_TRIANGLES)
        glVertex3f(10.0, 1.0, 1.0)
        glVertex3f(10.0, 1.0, -1.0)
        glVertex3f(13.0, 0.0, 0.0)

        glVertex3f(10.0, 1.0, -1.0)
        glVertex3f(10.0, -1.0, -1.0)
        glVertex3f(13.0, 0.0, 0.0)

        glVertex3f(10.0, -1.0, -1.0)
        glVertex3f(10.0, -1.0, 1.0)
        glVertex3f(13.0, 0.0, 0.0)

        glVertex3f(10.0, -1.0, 1.0)
        glVertex3f(10.0, 1.0, 1.0)
        glVertex3f(13.0, 0.0, 0.0)


        glVertex3f(1.0, 10.0, 1.0)
        glVertex3f(1.0, 10.0, -1.0)
        glVertex3f(0.0, 13.0, 0.0)

        glVertex3f(1.0, 10.0, -1.0)
        glVertex3f(-1.0, 10.0, -1.0)
        glVertex3f(0.0, 13.0, 0.0)

        glVertex3f(-1.0, 10.0, -1.0)
        glVertex3f(-1.0, 10.0, 1.0)
        glVertex3f(0.0, 13.0, 0.0)

        glVertex3f(-1.0, 10.0, 1.0)
        glVertex3f(1.0, 10.0, 1.0)
        glVertex3f(0.0, 13.0, 0.0)


        glVertex3f(1.0, 1.0, 10.0)
        glVertex3f(-1.0, 1.0, 10.0)
        glVertex3f(0.0, 0.0, 13.0)

        glVertex3f(-1.0, 1.0, 10.0)
        glVertex3f(-1.0, -1.0, 10.0)
        glVertex3f(0.0, 0.0, 13.0)

        glVertex3f(-1.0, -1.0, 10.0)
        glVertex3f(1.0, -1.0, 10.0)
        glVertex3f(0.0, 0.0, 13.0)

        glVertex3f(1.0, -1.0, 10.0)
        glVertex3f(1.0, 1.0, 10.0)
        glVertex3f(0.0, 0.0, 13.0)
        glEnd()

        # These loops draw each vertex on each face of each panel of
        # the detector.
        glBegin(GL_QUADS)
        glMaterial(GL_FRONT_AND_BACK, GL_AMBIENT, (1.0, 1.0, 0.8, 1.0))
        # For each panel
        for i in range(num_panels):
            # For each face excluding the one nearest the crystal.
            for j in [0, 2, 3, 4, 5]:
                glNormal3f(d_r_f[i][j][0, 0], d_r_f[i][j][0, 1],
                           d_r_f[i][j][0, 2])
                # For each vertex on the face
                for k in range(1, 5):
                    glVertex3f(d_r_f[i][j][k, 0], d_r_f[i][j][k, 1],
                               d_r_f[i][j][k, 2])
        glEnd()

        # Now drawing the front face of the panel, with the diffraction
        # image and goniometer shadow on it.
        glEnable(GL_TEXTURE_2D)

        ID = glGenTextures(1)
        glPixelStorei(GL_UNPACK_ALIGNMENT, 1)
        glBindTexture(GL_TEXTURE_2D, ID)
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT)
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT)
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR)
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR)
        glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL)

        if shadowsetting is True:
            glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, tex2.size[0], tex2.size[1],
            0, GL_LUMINANCE, GL_UNSIGNED_BYTE, tex_data)

        else:
            glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, tex2.size[0], tex2.size[1],
            0, GL_LUMINANCE, GL_UNSIGNED_BYTE, noshad_data)

        glBegin(GL_QUADS)
        glNormal3f(d_r_f[0][1][0, 0], d_r_f[0][1][0, 1], d_r_f[0][1][0, 2])
        glTexCoord2f(0.0, 0.0); glVertex3f(d_r_f[0][1][1, 0],
        d_r_f[0][1][1, 1], d_r_f[0][1][1, 2])
        glTexCoord2f(1.0, 0.0); glVertex3f(d_r_f[0][1][2, 0],
        d_r_f[0][1][2, 1], d_r_f[0][1][2, 2])
        glTexCoord2f(1.0, 1.0); glVertex3f(d_r_f[0][1][3, 0],
        d_r_f[0][1][3, 1], d_r_f[0][1][3, 2])
        glTexCoord2f(0.0, 1.0); glVertex3f(d_r_f[0][1][4, 0],
        d_r_f[0][1][4, 1], d_r_f[0][1][4, 2])
        glEnd()

        glDisable(GL_TEXTURE_2D)

        if self.size is None:
            self.size = self.GetClientSize()
        w, h = self.size
        w = max(w, 1.0)
        h = max(h, 1.0)
        xScale = 180.0 / w
        yScale = 180.0 / h
        glRotatef((self.y - self.lasty) * yScale, 1.0, 0.0, 0.0);
        glRotatef((self.x - self.lastx) * xScale, 0.0, 1.0, 0.0);

        self.SwapBuffers()

class MainWindow(wx.Frame):
    def __init__(self, parent = None, id = -1, title = "CBFviewer"):

        wx.Frame.__init__(
                self, parent, id, title, size = (1000, 1000),
                style = wx.DEFAULT_FRAME_STYLE | wx.NO_FULL_REPAINT_ON_RESIZE)

        box = wx.BoxSizer(wx.HORIZONTAL)
        box.Add(ViewerCanvas(self), 1, wx.EXPAND)

        self.SetAutoLayout(True)
        self.SetSizer(box)
        self.Layout()
        self.CreateStatusBar()

        # Creating the file menu
        filemenu = wx.Menu()
        menuitem = filemenu.Append(wx.ID_OPEN, "&Load", "Load a cbf image")
        self.Bind(wx.EVT_MENU, self.OnLoad, menuitem)
        menuitem = filemenu.Append(wx.ID_ABOUT, "&About",
                                   "Information about this program")
        self.Bind(wx.EVT_MENU, self.OnAbout, menuitem)
        filemenu.AppendSeparator()
        menuitem = filemenu.Append(wx.ID_EXIT, "E&xit",
                                   "Terminate the program")
        self.Bind(wx.EVT_MENU, self.OnExit, menuitem)

        # Creating the edit menu
        editmenu = wx.Menu()
        opsmenu = wx.Menu()
        menuitem = opsmenu.Append(-1, "Adjust &detector\tCtrl+D",
                                  "Bring up detector orientation dialog")
        self.Bind(wx.EVT_MENU, self.OnAdjustDet, menuitem)
        menuitem = opsmenu.Append(-1, "Adjust &sample\tCtrl+S",
                                  "Bring up sample orientation dialog")
        self.Bind(wx.EVT_MENU, self.OnAdjustSam, menuitem)
        menuitem = opsmenu.Append(-1, "Adjust &goniometer\tCtrl+G",
                                  "Bring up goniometer orientation dialog")
        self.Bind(wx.EVT_MENU, self.OnAdjustGonio, menuitem)
        editmenu.AppendMenu(wx.ID_PREFERENCES, "&Options", opsmenu)
        menuitem = editmenu.Append(-1, "&Reset\tCtrl+R",
                                   "Reset parameters to original values")
        self.Bind(wx.EVT_MENU, self.OnReset, menuitem)
        menuitem = editmenu.Append(-1, "&Shadow graph",
                                   "Bring up affected pixel plotting dialog")
        self.Bind(wx.EVT_MENU, self.OnPixPlot, menuitem)

        # Creating the view menu
        viewmenu = wx.Menu()
        menuitem = viewmenu.Append(-1, "x-y\tCtrl+1",
                    "Change viewing projection to be parallel to x-y plane")
        self.Bind(wx.EVT_MENU, self.xy, menuitem)
        menuitem = viewmenu.Append(-1, "x-z\tCtrl+2",
                    "Change viewing projection to be parallel to x-z plane")
        self.Bind(wx.EVT_MENU, self.xz, menuitem)
        menuitem = viewmenu.Append(-1, "y-z\tCtrl+3",
                    "Change viewing projection to be parallel to y-z plane")
        self.Bind(wx.EVT_MENU, self.yz, menuitem)
        menuitem = viewmenu.Append(-1, "y-x\tCtrl+4",
                    "Change viewing projection to be parallel to y-x plane")
        self.Bind(wx.EVT_MENU, self.yx, menuitem)
        menuitem = viewmenu.Append(-1, "z-x\tCtrl+5",
                    "Change viewing projection to be parallel to z-x plane")
        self.Bind(wx.EVT_MENU, self.zx, menuitem)
        menuitem = viewmenu.Append(-1, "z-y\tCtrl+6",
                    "Change viewing projection to be parallel to z-y plane")
        self.Bind(wx.EVT_MENU, self.zy, menuitem)

        # Creating the menubar
        menubar = wx.MenuBar()
        menubar.Append(filemenu,"&File")
        menubar.Append(editmenu, "&Edit")
        menubar.Append(viewmenu, "&View")
        self.SetMenuBar(menubar)

        self.Show(True)
        self.Centre()

    def xy(self, event):
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()
        glRotatef(180.0, 0.0, 1.0, 0.0)

        self.Refresh(False)

    def xz(self, event):
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()
        glRotatef(90.0, 1.0, 0.0, 0.0)
        glRotatef(180.0, 0.0, 0.0, 1.0)

        self.Refresh(False)

    def yz(self, event):
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()
        glRotatef(90.0, 1.0, 0.0, 0.0)
        glRotatef(-90.0, 0.0, 0.0, 1.0)

        self.Refresh(False)

    def yx(self, event):
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()
        glRotatef(-90.0, 0.0, 0.0, 1.0)

        self.Refresh(False)

    def zx(self, event):
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()
        glRotatef(-90.0, 0.0, 1.0, 0.0)
        glRotatef(90.0, 0.0, 0.0, 1.0)

        self.Refresh(False)

    def zy(self, event):
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()
        glRotatef(-90.0, 0.0, 1.0, 0.0)

        self.Refresh(False)

    def OnLoad(self, event):
        wildcard = "cbf images (*.cbf)|*.cbf|All Files (*.*)|*.*"
        dlg = wx.FileDialog(self, message = "Select a cbf image",
                            defaultDir = os.getcwd(), wildcard = wildcard,
                            defaultFile = "", style = wx.OPEN)
        done = False
        while not done:
            try:
                if dlg.ShowModal() == wx.ID_OK:
                    filename = dlg.GetPath()
                    self.Close(True)
                    launcher.view_experiment(filename)
                    done = True
                elif dlg.ShowModal() == wx.ID_CANCEL:
                    done = True
                    dlg.Destroy()

            except IOError:
                message = "Invalid file. Please choose another."
                caption = "Error"
                wx.MessageBox(message, caption, wx.ICON_ERROR)


    def OnAdjustDet(self, event):

        det_angle = DetectorDialog(None, title='Detector Orientation')
        det_angle.ShowModal()

    def OnAdjustSam(self, event):

        sam_angles = SampleDialog(None, title='Sample Orientation')
        sam_angles.ShowModal()

    def OnAdjustGonio(self, event):

        gonio_angles = GonioDialog(None, title='Goniometer Orientation')
        gonio_angles.ShowModal()

    def OnReset(self, event):
        self.Close(True)
        launcher.view_experiment(filename)

    def OnPixPlot(self, event):
        pix_plotter = PlotDialog(None, title='Graph Plotting')
        pix_plotter.ShowModal()

    def OnAbout(self,event):
        message = ("CBFViewer uses PyOpenGL in wxPython to display a"
                   "3D model of diffraction experiments")
        caption = "About CBFviewer"
        wx.MessageBox(message, caption, wx.OK)

    def OnExit(self,event):
        self.Close(True)

class DetectorDialog(wx.Dialog):

    def __init__(self, *args, **kw):
        super(DetectorDialog, self).__init__(*args, **kw)
        self.InitUI()
        self.SetSize((350, 310))
        self.SetTitle("Adjust detector orientation")
        self.Centre()

    def InitUI(self):

        panel = wx.Panel(self)
        wx.StaticText(panel, label = "Enter a new detector angle and "
                      "distance.\nAngles are measured in degrees "
                      "anti-\nclockwise from the positive z axis, and "
                      "\ndistances in milimetres from the sample.",
                      pos = (20,20))
        wx.StaticText(panel, label = 'New angle: ', pos = (20, 100))
        wx.StaticText(panel, label = 'New distance: ', pos = (20, 150))

        # angles widget
        self.angles = fs.FloatSpin(panel, -1, pos=(150, 95), min_val = -90.0,
                                   max_val = 90.0, increment = 0.5,
                                   value = round(DET_ANGLE, 1),
                                   agwStyle = fs.FS_LEFT)
        self.angles.SetFormat("%f")
        self.angles.SetDigits(1)
        # distance widget
        self.dist = fs.FloatSpin(panel, -1, pos=(150, 145), min_val = 20.0,
                                 max_val = 500.0, increment = 1.0,
                                 value = round(d, 1),
                                 agwStyle = fs.FS_LEFT)
        self.dist.SetFormat("%f")
        self.dist.SetDigits(1)
        # detector texture toggle
        self.detcheck = wx.CheckBox(panel, -1, 'Enable detector texture',
                                    (40, 200))
        self.detcheck.SetValue((not blanksetting))

        #obtn = wx.Button(panel, label='&Ok', pos=(30, 230))
        rbtn = wx.Button(panel, label='&Adjust', pos=(70, 230))
        rbtn.SetFocus()
        cbtn = wx.Button(panel, wx.ID_CLOSE, label='&Close', pos=(185, 230))

        #obtn.Bind(wx.EVT_BUTTON, self.OnOk)
        rbtn.Bind(wx.EVT_BUTTON, self.OnAdjust)
        cbtn.Bind(wx.EVT_BUTTON, self.OnClose)

    def OnAdjust(self, event):

        global d
        global DET_ANGLE
        global tex_data
        global tex2
        global tex_noshad
        global noshad_data
        global blanksetting
        global tex
        global tex_shadow
        # angle = angular displacement
        angle = np.deg2rad(self.angles.GetValue() - DET_ANGLE)
        displacement = self.dist.GetValue() - d
        # Rotation matrix
        r_x = np.array([[1.0, 0.0, 0.0],
                        [0.0, np.cos(angle), -np.sin(angle)],
                        [0.0, np.sin(angle), np.cos(angle)]], float)
        det = []

        for i in range(num_panels):
            det.append(modeller.Panel(dimensions[i, 0], dimensions[i, 1]))
            detector_origin[i] = (detector_origin[i] + displacement
                                  * detector_normal)
            detector_origin[i] = np.dot(r_x, detector_origin[i])
            fast_axis[i] = np.dot(r_x, fast_axis[i])
            slow_axis[i] = np.dot(r_x, slow_axis[i])
            d_r_f[i] = det[i].get_coords(detector_origin[i], fast_axis[i],
                                         slow_axis[i])

        d = self.dist.GetValue()
        DET_ANGLE = self.angles.GetValue()
        blanksetting = not self.detcheck.GetValue()

        newshadow = shadowmapper.Shadow(shadowmapper.day2_coords, filename,
        DET_ANGLE, d, detector_origin, fast_axis, slow_axis)

        [tex_data, tex, tex2, tex_shadow, tex_noshad,
        noshad_data] = newshadow.shadow_texture(kappa = kappa, omega = omega,
                                                blank = blanksetting)

        pixels = newshadow.affected_pixels(kappa = kappa, omega = omega,
                                          blank = blanksetting)[1]

        if shadowsetting is True:
            print 'At these goniometer angles, the shadow affects', pixels, \
                  'pixels.'

        self.Destroy()
        app.MainLoop()

    #def OnOk(self, event):
        #self.OnApply(evt)
        #self.Destroy()

    def OnClose(self, event):
        self.Destroy()

class SampleDialog(wx.Dialog):

    def __init__(self, *args, **kw):
        super(SampleDialog, self).__init__(*args, **kw)

        self.InitUI()
        self.SetSize((350, 380))
        self.SetTitle("Adjust sample orientation")
        self.Centre()

    def InitUI(self):

        self.angle = dict()
        self.rb = dict()
        panel = wx.Panel(self)
        wx.StaticText(panel, label = "Enter new sample angles between "
                      "-90.0 and 90.0 degrees. ", pos = (20,20))

        labels = ['Rotate sample about x axis (omega): ',
        'Rotate sample about y axis: ', 'Rotate sample about z axis: ',
        'Sample position: ']

        for i in range(4):
            wx.StaticText(panel, label = labels[i], pos = (20, (80 + 50 * i)))

        for i in range(3):          # creating 3 FloatSpin widgets
            self.angle[i] = fs.FloatSpin(panel, -1, pos=(220, (75 + 50 * i)),
                                         min_val = -90.0, max_val = 90.0,
                                         increment = 0.5,
                                         value = round(sample_angles[i], 1),
                                         agwStyle = fs.FS_LEFT)
            self.angle[i].SetFormat("%f")
            self.angle[i].SetDigits(1)

        self.rb1 = wx.RadioButton(panel, -1, '(0, 0, 0)', (220, 210),
                                  style = wx.RB_GROUP)
        self.rb2 = wx.RadioButton(panel, -1, '(6.4, 0, 0)', (220, 230))
        self.rb3 = wx.RadioButton(panel, -1, '(-6.4, 0, 0)', (220, 250))

        rbtn = wx.Button(panel, label='&Adjust', pos=(70, 300))
        rbtn.SetFocus()
        cbtn = wx.Button(panel, wx.ID_CLOSE, label='&Close', pos=(185, 300))

        rbtn.Bind(wx.EVT_BUTTON, self.OnCompute)
        cbtn.Bind(wx.EVT_BUTTON, self.OnClose)

    def OnCompute(self, event):

        alpha = -np.deg2rad(self.angle[0].GetValue())
        beta  = np.deg2rad(self.angle[1].GetValue())
        gamma = np.deg2rad(self.angle[2].GetValue())

        # rotation matrices for each axis
        r_x = np.array([[1.0, 0.0, 0.0],
        [0.0, np.cos(alpha), -np.sin(alpha)],
        [0.0, np.sin(alpha), np.cos(alpha)]], float)
        r_y = np.array([[np.cos(beta), 0.0, np.sin(beta)],
        [0.0, 1.0, 0.0], [-np.sin(beta), 0.0, np.cos(beta)]], float)
        r_z = np.array([[np.cos(gamma), -np.sin(gamma), 0.0],
        [np.sin(gamma), np.cos(gamma), 0.0], [0.0, 0.0, 1.0]], float)
        # r_tot = matrix representing the total rotation
        r_tot = np.dot(r_z, np.dot(r_y, r_x))

        newSample = modeller.Sample()
        r_f = newSample.rotate(r_tot[0], r_tot[1], r_tot[2])

        sample_x = 0.0
        if self.rb2.GetValue() is True:
            sample_x = 6.4
        elif self.rb3.GetValue() is True:
            sample_x = -6.4

        new_origin = deepcopy(detector_origin)
        new_origin[0, 0] += sample_x

        newDet = []
        for i in range(num_panels):
            newDet.append(modeller.Panel(dimensions[i, 0], dimensions[i, 1]))
            d_r_f[i] = (newDet[i].get_coords(new_origin, fast_axis[i],
                        slow_axis[i]))

        print 'The new alpha, beta, and gamma angles are, respectively, ', \
               np.rad2deg(alpha), ', ', np.rad2deg(beta), 'and ', \
               np.rad2deg(gamma), 'degrees'

        newshadow = shadowmapper.Shadow(shadowmapper.day2_coords, filename,
        DET_ANGLE, d, detector_origin, fast_axis, slow_axis)

        [tex_data, tex, tex2, tex_shadow, tex_noshad,
        noshad_data] = newshadow.shadow_texture(kappa = kappa, omega = omega, blank = blanksetting)

        pixels = newshadow.affected_pixels(kappa = kappa, omega = omega, blank = blanksetting)[1]

        if shadowsetting is True:
            print 'At these goniometer angles, the shadow affects', pixels, \
                  'pixels.'

        self.Destroy()
        app.MainLoop()

    def OnClose(self, event):

        self.Destroy()

class GonioDialog(wx.Dialog):

    def __init__(self, *args, **kw):
        super(GonioDialog, self).__init__(*args, **kw)

        self.InitUI()
        self.SetSize((350, 500))
        self.SetTitle("Adjust goniometer angles")
        self.Centre()

    def InitUI(self):

        self.widget = dict()
        panel = wx.Panel(self)
        wx.StaticText(panel, label = 'Enter new goniometer parameters.',
                      pos = (20,20))

        labels = ['Adjust omega: ', 'Adjust kappa: ', 'x offset: ',
                  'y offset: ', 'z offset: ']
        values = [round(omega, 1), round(kappa, 1), offset_matrix[0],
                  offset_matrix[1], offset_matrix[2]]
        minmaxinc = np.array([[-270.0, 270.0, 0.5], [-190.0, 10.0, 0.5],
                  [-2.5, 2.5, 0.1], [-2.5, 2.5, 0.1], [-2.5, 2.5, 0.1]])

        # Creating 5 FloatSpin widgets and their labels.
        for i in range(5):
            wx.StaticText(panel, label = labels[i], pos = (20, (80 + 50 * i)))
            self.widget[i] = fs.FloatSpin(panel, -1, pos=(200, (75 + 50 * i)),
            min_val = minmaxinc[i, 0], max_val = minmaxinc[i, 1],
            increment = minmaxinc[i, 2], value = round(values[i], 1),
            agwStyle = fs.FS_LEFT)
            self.widget[i].SetFormat("%f")
            self.widget[i].SetDigits(1)

        self.shadowcheck = wx.CheckBox(panel, -1, 'Enable goniometer shadow',
                                      (40, 350))
        self.shadowcheck.SetValue(shadowsetting)

        rbtn = wx.Button(panel, label='&Adjust', pos=(70, 400))
        rbtn.SetFocus()
        cbtn = wx.Button(panel, wx.ID_CLOSE, label='&Close', pos=(185, 400))

        rbtn.Bind(wx.EVT_BUTTON, self.OnCompute)
        cbtn.Bind(wx.EVT_BUTTON, self.OnClose)

    def OnCompute(self, event):

        global tex_data
        global tex2
        global omega
        global kappa
        global offset_matrix
        global shadowsetting
        global tex_shadow
        global tex_noshad
        global noshad_data
        global tex

        omega = self.widget[0].GetValue()
        kappa = self.widget[1].GetValue()
        offset_matrix = [self.widget[2].GetValue(), self.widget[3].GetValue(),
                         self.widget[4].GetValue()]

        newshadow = shadowmapper.Shadow(shadowmapper.day2_coords, filename,
        DET_ANGLE, d, detector_origin, fast_axis, slow_axis)

        [tex_data, tex, tex2, tex_shadow, tex_noshad,
        noshad_data] = newshadow.shadow_texture(kappa = kappa, omega = omega, blank = blanksetting)

        pixels = newshadow.affected_pixels(kappa = kappa, omega = omega, blank = blanksetting)[1]

        shadowsetting = self.shadowcheck.GetValue()

        if shadowsetting is True:
            print 'At these goniometer angles, the shadow affects', pixels, \
                  'pixels.'

        self.Destroy()
        app.MainLoop()

    def OnClose(self, event):

        self.Destroy()

class PlotDialog(wx.Dialog):

    def __init__(self, *args, **kw):
        super(PlotDialog, self).__init__(*args, **kw)
        self.InitUI()
        self.SetSize((350, 170))
        self.SetTitle("Plot graph of affected pixels")
        self.Centre()

    def InitUI(self):

        panel = wx.Panel(self)
        wx.StaticText(panel, label = "Enter the number of points in "
                                      "each dimension. \n Note the "
                                      "plotting will take some time.",
                                      pos = (20,20))
        wx.StaticText(panel, label = "Number of points: ", pos = (20, 80))

        # angles widget
        self.points = wx.SpinCtrl(panel, -1, pos=(150, 75), min = 10,
                                  max = 10000, initial = 50)

        rbtn = wx.Button(panel, label='&Plot', pos=(70, 130))
        rbtn.SetFocus()
        cbtn = wx.Button(panel, wx.ID_CLOSE, label='&Close', pos=(185, 130))

        #obtn.Bind(wx.EVT_BUTTON, self.OnOk)
        rbtn.Bind(wx.EVT_BUTTON, self.OnPlot)
        cbtn.Bind(wx.EVT_BUTTON, self.OnClose)

    def OnPlot(self, event):
        myshadow = shadowmapper.Shadow(shadowmapper.day2_coords, filename,
                                      DET_ANGLE, d, detector_origin,
                                      fast_axis, slow_axis)
        myshadow.pixel_plotter(num_points = self.points.GetValue())
        self.Destroy()

    def OnClose(self, event):
        self.Destroy()

if __name__ == "__main__":
      # This is old test code and probably doesn't work anymore.
      detector_normal = (0.0, -0.4226182617406947, -0.9063077870366522)
      crystal_i = [0.7065882430629097, 0.5534899064860146, -0.440887716072214]
      crystal_j = [0.4514320564090011, -0.832380862760658, -0.321482810980872]
      crystal_k = [-0.5449239884714299, 0.02812512627395, -0.8380136180638496]

      detector = modeller.Detector()

      d_r_f = detector.get_coords(detector_normal)

      crystal = modeller.Sample()

      r_f = crystal.rotate(crystal_i, crystal_j, crystal_k)

      B = (0.0, 0.0, -1.0)

      app = wx.App()
      frame = MainWindow()
      app.MainLoop()
