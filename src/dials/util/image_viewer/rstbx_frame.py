from __future__ import annotations

import os
import pickle

import wx
import wx.lib.colourselect

import rstbx.viewer.display
import wxtbx.plots
from libtbx.utils import to_unicode
from wxtbx import bitmaps, icons

from dials.util import Sorry

# Instance to bind external update event to an event handler
EVT_EXTERNAL_UPDATE = wx.PyEventBinder(wx.NewEventType(), 0)


class ExternalUpdateEvent(wx.PyCommandEvent):
    """XXX This class, along with the EVT_EXTERNAL_UPDATE instance
    should perhaps move into its own file?
    """

    def __init__(self, eventType=EVT_EXTERNAL_UPDATE.evtType[0], id=0):
        wx.PyCommandEvent.__init__(self, eventType, id)
        self.img = None
        self.title = None


class XrayFrame(wx.Frame):
    # Maximum number of entries in the chooser.
    CHOOSER_SIZE = 1024

    def __init__(self, *args, **kwds):
        super().__init__(*args, **kwds)
        self.settings = rstbx.viewer.settings()
        self.viewer = rstbx.viewer.display.XrayView(self, -1, size=(1024, 640))
        self.viewer.SetMinSize((640, 640))
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(self.sizer)
        self.sizer.Add(self.viewer, 1, wx.EXPAND)
        self.statusbar = self.CreateStatusBar()
        self.settings_frame = None
        self.zoom_frame = None
        self.zoom_3d = None
        self.plot_frame = None
        self._img = None
        self._distl = None
        self.toolbar = self.CreateToolBar(style=wx.TB_TEXT)
        self.setup_toolbar()
        self.toolbar.Realize()
        self.mb = wx.MenuBar()
        self.setup_menus()
        self.SetMenuBar(self.mb)
        self.Fit()
        self.SetMinSize(self.GetSize())
        self.OnShowSettings(None)
        self.Bind(EVT_EXTERNAL_UPDATE, self.OnExternalUpdate)

    def OnExternalUpdate(self, event):
        """The OnExternalUpdate() function updates the image and the title
        from @p event.
        """

        # See self.load_image().
        self._img = event.img
        self.viewer.set_image(self._img)
        if self.settings_frame is not None:
            self.settings_frame.set_image(self._img)
        if self.zoom_frame is not None:
            self.zoom_frame.set_image(self._img)
        self.SetTitle(event.title)
        self.update_statusbar()
        self.Layout()

    def setup_toolbar(self):
        btn = self.toolbar.AddTool(
            toolId=-1,
            label="Load file",
            bitmap=icons.hkl_file.GetBitmap(),
            shortHelp="Load file",
            kind=wx.ITEM_NORMAL,
        )
        self.Bind(wx.EVT_MENU, self.OnLoadFile, btn)
        btn = self.toolbar.AddTool(
            toolId=-1,
            label="Settings",
            bitmap=icons.advancedsettings.GetBitmap(),
            shortHelp="Settings",
            kind=wx.ITEM_NORMAL,
        )
        self.Bind(wx.EVT_MENU, self.OnShowSettings, btn)
        btn = self.toolbar.AddTool(
            toolId=-1,
            label="Zoom",
            bitmap=icons.search.GetBitmap(),
            shortHelp="Zoom",
            kind=wx.ITEM_NORMAL,
        )
        self.Bind(wx.EVT_MENU, self.OnZoom, btn)
        txt = wx.StaticText(self.toolbar, -1, "Image:")
        self.toolbar.AddControl(txt)
        self.image_chooser = wx.Choice(self.toolbar, -1, size=(300, -1))
        self.toolbar.AddControl(self.image_chooser)
        self.Bind(wx.EVT_CHOICE, self.OnChooseImage, self.image_chooser)
        btn = self.toolbar.AddTool(
            toolId=wx.ID_BACKWARD,
            label="Previous",
            bitmap=bitmaps.fetch_icon_bitmap("actions", "1leftarrow"),
            shortHelp="Previous",
            kind=wx.ITEM_NORMAL,
        )
        self.Bind(wx.EVT_MENU, self.OnPrevious, btn)
        btn = self.toolbar.AddTool(
            toolId=wx.ID_FORWARD,
            label="Next",
            bitmap=bitmaps.fetch_icon_bitmap("actions", "1rightarrow"),
            shortHelp="Next",
            kind=wx.ITEM_NORMAL,
        )
        self.Bind(wx.EVT_MENU, self.OnNext, btn)

    def setup_menus(self):
        file_menu = wx.Menu()
        self.mb.Append(file_menu, "File")
        item = file_menu.Append(-1, "Open integration results...")
        self.Bind(wx.EVT_MENU, self.OnLoadIntegration, item)
        item = file_menu.Append(-1, "Open image...")
        self.Bind(wx.EVT_MENU, self.OnLoadFile, item)
        actions_menu = wx.Menu()
        self.mb.Append(actions_menu, "Actions")
        item = actions_menu.Append(-1, "Change beam center...")
        self.Bind(wx.EVT_MENU, self.OnChangeBeamCenter, item)
        item = actions_menu.Append(-1, "Reset beam center to header value")
        self.Bind(wx.EVT_MENU, lambda evt: self.viewer.ResetBeamCenter(), item)
        item = actions_menu.Append(-1, "Save screenshot...")
        self.Bind(wx.EVT_MENU, self.OnScreenShot, item)

    def get_key(self, file_name_or_data):
        """The get_key() function returns the key of @p file_name_or_data.
        In the case of dictionaries, it is the timestamp of the image.
        For file names, the key is an ASCII-encoded absolute path string.
        Otherwise, get_key() returns @c None.
        """

        try:
            return file_name_or_data["TIMESTAMP"]
        except TypeError:
            pass

        try:
            return file_name_or_data.get_image_file()
        except AttributeError:
            pass

        try:
            return os.path.abspath(file_name_or_data)
        except TypeError:
            pass

        try:
            return os.path.abspath(file_name_or_data.encode("ascii"))
        except (AttributeError, TypeError):
            pass

        return None

    def load_image(self, file_name_or_data):
        """The load_image() function displays the image from @p
        file_name_or_data.  The chooser is updated appropriately.
        """

        key = self.get_key(file_name_or_data)
        if isinstance(file_name_or_data, dict):
            self._img = rstbx.viewer.image(file_name_or_data)
        else:
            try:
                self._img = rstbx.viewer.image(key)
            except OSError:
                raise Sorry(
                    (
                        "The file '%s' could not be recognized as a supported "
                        + "image format; please make sure it is actually a detector image."
                    )
                    % key
                )

        # Update the selection in the chooser.
        i = self.add_file_name_or_data(file_name_or_data)
        self.image_chooser.SetSelection(i)

        self.viewer.set_image(self._img)
        if self.settings_frame is not None:
            self.settings_frame.set_image(self._img)
        self.SetTitle(to_unicode(key))
        self.update_statusbar()
        self.Layout()

    def load_distl_output(self, file_name):
        with open(file_name, "rb") as fh:
            distl = pickle.load(fh)
        self._distl = distl
        img_files = []
        for img_id in sorted(distl.images.keys()):
            img = distl.images[img_id]
            img_files.append(img["relpath"])
        if len(img_files) == 0:
            raise Sorry("No images in this result!")
        self.image_chooser.SetItems([os.path.basename(f) for f in img_files])
        self.image_chooser.SetSelection(0)
        self.load_image(img_files[0])
        self.annotate_image(img_files[0])

    def add_file_name_or_data(self, file_name_or_data):
        """The add_file_name_or_data() function appends @p
        file_name_or_data to the image chooser, unless it is already
        present.  For file-backed images, the base name is displayed in
        the chooser.  If necessary, the number of entries in the chooser
        is pruned.  The function returns the index of the recently added
        entry.  XXX This is probably the place for heuristics to determine
        if the viewer was given a pattern, or a plain list of files.  XXX
        Rename this function, because it only deals with the chooser?
        """

        key = self.get_key(file_name_or_data)
        for i in range(self.image_chooser.GetCount()):
            if key == self.image_chooser.GetClientData(i):
                return i
        if self.image_chooser.GetCount() >= self.CHOOSER_SIZE:
            self.image_chooser.Delete(0)
        i = self.image_chooser.GetCount()
        if isinstance(file_name_or_data, dict):
            self.image_chooser.Insert(key, i, None)
        else:
            self.image_chooser.Insert(os.path.basename(key), i, key)
        return i

    def annotate_image(self, file_name):
        assert self._distl is not None
        for img_id in sorted(self._distl.images.keys()):
            img = self._distl.images[img_id]
            if img["relpath"] == file_name:
                spots = img["spotoutput"]["inlier_spots"]
                self._img.set_spots(spots)
                break

    def load_integration(self, dir_name):
        from rstbx.viewer import processing

        assert os.path.isdir(dir_name)
        self.proc_frame = processing.ProcessingFrame(None, -1, "LABELIT")
        self.proc_frame.set_viewer_frame(self)
        self.proc_frame.LoadResults(dir_name)
        self.proc_frame.Show()

    def display_integration_result(self, result):
        assert isinstance(result, dict)
        self._img.set_integration_results(result)
        self.viewer.Refresh()

    def update_statusbar(self, info=None):
        if info is None:
            self.statusbar.SetStatusText(
                "Click and drag to pan, mouse wheel or double click to zoom"
            )
        else:
            self.statusbar.SetStatusText(info.format())

    def update_settings(self, layout=True):
        self.viewer.update_settings(layout)

    def set_brightness(self, brightness):
        if (brightness > 0) and (brightness <= 500):
            self.settings.brightness = brightness
            if self.settings_frame is not None:
                self.settings_frame.update_controls()
            self.viewer.update_settings(layout=False)

    def OnLoadFile(self, event):
        wildcard_str = ""
        if wx.PlatformInfo[4] != "wxOSX-cocoa":
            from iotbx import file_reader

            wildcard_str = file_reader.get_wildcard_string("img")
        file_name = wx.FileSelector(
            "Image file",
            wildcard=wildcard_str,
            default_path="",
            flags=wx.FD_OPEN,
        )
        if file_name != "":
            self.load_image(file_name)

    def OnLoadLabelitResult(self, event):
        file_name = wx.FileSelector("Labelit result", default_path="", flags=wx.FD_OPEN)
        if file_name != "":
            self.load_image(file_name)

    def OnLoadIntegration(self, event):
        dir_name = wx.DirSelector("Integration result", defaultPath="")
        if dir_name != "":
            self.load_integration(dir_name)

    def OnShowSettings(self, event):
        raise NotImplementedError("Removed due to non-used code path")

    def OnShowZoom(self, event):
        if self.zoom_frame is None:
            self.zoom_frame = ZoomFrame(
                self,
                -1,
                "Zoom",
                style=wx.CAPTION | wx.CLOSE_BOX | wx.RESIZE_BORDER | wx.SYSTEM_MENU,
            )
            self.zoom_frame.set_image(self._img)
            self.zoom_frame.Show()
        self.zoom_frame.Raise()

    def OnShow3D(self, event):
        if self.zoom_3d is None:
            from rstbx.viewer import pixels3d

            self.zoom_3d = pixels3d.pixel_viewer_3d_frame(
                self,
                -1,
                "3D view",
                style=wx.CAPTION | wx.CLOSE_BOX | wx.RESIZE_BORDER | wx.SYSTEM_MENU,
            )
            self.zoom_3d.set_image(self._img)
            self.zoom_3d.Show()
        self.zoom_3d.Raise()

    def OnShowPlot(self, event):
        if self.plot_frame is None:
            self.plot_frame = PlotFrame(
                self,
                -1,
                "Intensity profile",
                style=wx.CAPTION | wx.CLOSE_BOX | wx.SYSTEM_MENU,
            )
            self.plot_frame.Show()

    def OnZoom(self, event):
        if self.settings.zoom_level == 6:
            self.settings.zoom_level = 0
        else:
            self.settings.zoom_level += 1
        if hasattr(self.viewer, "update_settings"):
            self.viewer.update_settings(layout=True)
        if self.settings_frame is not None:
            self.settings_frame.update_controls()

    def OnChooseImage(self, event):
        self.load_image(
            self.image_chooser.GetClientData(self.image_chooser.GetSelection())
        )

    def OnPrevious(self, event):
        n = self.image_chooser.GetSelection()
        if n != wx.NOT_FOUND and n - 1 >= 0:
            self.load_image(self.image_chooser.GetClientData(n - 1))

    def OnNext(self, event):
        n = self.image_chooser.GetSelection()
        if n != wx.NOT_FOUND and n + 1 < self.image_chooser.GetCount():
            self.load_image(self.image_chooser.GetClientData(n + 1))

    def OnScreenShot(self, event):
        file_name = wx.FileSelector(
            default_filename="xray.png",
            default_path="",
            wildcard="PNG image (*.png)|*.png",
            flags=wx.SAVE,
        )
        if file_name != "":
            self.viewer.save_image(file_name)

    def OnChangeBeamCenter(self, event):
        wx.MessageBox("Click on any point in the image to set the new beam center.")
        self.statusbar.SetStatusText("Changing beam center")
        self.viewer.ChangeBeamCenter()


mag_levels = [8, 16, 24, 32, 48, 64]


class ZoomFrame(wx.MiniFrame):
    def __init__(self, *args, **kwds):
        super().__init__(*args, **kwds)
        self.settings = self.GetParent().settings
        self.control_panel = wx.Panel(self)
        self.panel = rstbx.viewer.display.ZoomView(self, -1)
        szr = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(szr)
        szr.Add(self.control_panel)
        szr.Add(self.panel, 1, wx.EXPAND)
        self.numbers_box = wx.CheckBox(self.control_panel, -1, "Show intensity values")
        txt1 = wx.StaticText(self.control_panel, -1, "Text color:")
        self.text_color = wx.lib.colourselect.ColourSelect(
            self.control_panel, colour=(255, 255, 0)
        )
        pszr = wx.BoxSizer(wx.VERTICAL)
        self.control_panel.SetSizer(pszr)
        box1 = wx.BoxSizer(wx.HORIZONTAL)
        pszr.Add(box1)
        box1.Add(self.numbers_box, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        box1.Add(txt1, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        box1.Add(self.text_color, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        box2 = wx.BoxSizer(wx.HORIZONTAL)
        pszr.Add(box2)
        txt2 = wx.StaticText(self.control_panel, -1, "Magnification:")
        self.mag_ctrl = wx.Choice(
            self.control_panel, -1, choices=["%dx" % x for x in mag_levels]
        )
        self.mag_ctrl.SetSelection(1)
        box2.Add(txt2, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        box2.Add(self.mag_ctrl, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        self.Bind(wx.EVT_CLOSE, lambda evt: self.Destroy(), self)
        self.Bind(wx.EVT_WINDOW_DESTROY, self.OnDestroy)
        self.Bind(wx.EVT_CHECKBOX, self.OnChangeSettings, self.numbers_box)
        self.Bind(
            wx.lib.colourselect.EVT_COLOURSELECT, self.OnChangeSettings, self.text_color
        )
        self.Bind(wx.EVT_CHOICE, self.OnChangeSettings, self.mag_ctrl)
        szr.Fit(self.panel)
        self.Fit()

    def __getattr__(self, name):
        return getattr(self.panel, name)

    def OnDestroy(self, event):
        self.GetParent().zoom_frame = None

    def OnChangeSettings(self, event):
        self.panel.flag_show_intensities = self.numbers_box.GetValue()
        self.panel.text_color = self.text_color.GetValue()
        zoom = mag_levels[self.mag_ctrl.GetSelection()]
        self.panel.zoom_level = zoom
        self.Refresh()


class PlotFrame(wx.MiniFrame):
    def __init__(self, *args, **kwds):
        super().__init__(*args, **kwds)
        self.plot = LinePlot(self, figure_size=(8, 3))
        szr = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(szr)
        szr.Add(self.plot, 1, wx.EXPAND)
        self.Fit()
        self.Bind(wx.EVT_CLOSE, lambda evt: self.Destroy(), self)
        self.Bind(wx.EVT_WINDOW_DESTROY, self.OnDestroy)

    def __getattr__(self, name):
        return getattr(self.plot, name)

    def OnDestroy(self, event):
        self.GetParent().plot_frame = None


class LinePlot(wxtbx.plots.plot_container):
    def show_plot(self, line):
        self.figure.clear()
        ax = self.figure.add_subplot(111)
        x_data = list(range(len(line.values)))
        ax.plot(x_data, line.values, "b-", linewidth=1)
        ax.set_ylabel("Intensity")
        if line.lattice_length is not None:
            ax.set_title(
                f"Line distance = {line.distance:.2f}mm; avg. lattice length = {line.lattice_length:.2f} Angstrom"
            )
        else:
            ax.set_title(f"Line distance: {line.distance:.2f}mm")
        self.canvas.draw()
