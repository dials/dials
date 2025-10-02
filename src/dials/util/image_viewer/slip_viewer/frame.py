#
# Known issues: Recentering on resize and when switching between
# different image types.  Ring centre on image switch.


from __future__ import annotations

import math
import os
from importlib.machinery import SourceFileLoader

import wx

from boost_adaptbx.boost.python import c_sizeof
from rstbx.viewer import image as rv_image
from rstbx.viewer import settings as rv_settings
from wxtbx import bitmaps

from ..rstbx_frame import EVT_EXTERNAL_UPDATE
from ..rstbx_frame import XrayFrame as XFBaseClass
from . import pyslip, tile_generation
from .calibration_frame import SBSettingsFrame

pyslip._Tiles = tile_generation._Tiles

# Use minimum value for an int to indicate a masked pixel
int_bits = c_sizeof("int") * 8
MASK_VAL = -(2 ** (int_bits - 1))


class chooser_wrapper:
    def __init__(self, image_set, index):
        self.image_set = image_set
        self.path = os.path.basename(image_set.get_path(index))
        self.full_path = image_set.get_path(index)
        self.index = index
        self._image_data = None

    def __str__(self):
        return "%s [%d]" % (self.path, self.index + 1)

    def get_detector(self):
        return self.image_set.get_detector()

    def get_scan(self):
        return self.image_set.get_scan()

    def get_beam(self):
        return self.image_set.get_beam()

    def get_mask(self):
        return self.image_set.get_mask(self.index)

    def get_image_data(self, corrected=True):
        if self._image_data is None:
            if corrected:
                return self.image_set.get_corrected_data(self.index)
            else:
                return self.image_set.get_raw_data(self.index)
        return self._image_data

    def set_image_data(self, image_data):
        self._image_data = image_data

    def get_detectorbase(self):
        return self.image_set.get_detectorbase(self.index)

    def get_vendortype(self):
        return self.image_set.get_vendortype(self.index)

    def show_header(self):
        return self.image_set.get_detectorbase(self.index).show_header()


class XrayFrame(XFBaseClass):
    def set_pyslip(self, parent):
        self.pyslip = pyslip.PySlip(parent, tile_dir=None, min_level=0)

    def __init__(self, *args, **kwds):
        self.params = kwds.get("params", None)
        if "params" in kwds:
            del kwds["params"]  # otherwise wx complains

        ### Collect any plugins
        slip_viewer_dir = os.path.join(os.path.dirname(__file__))
        contents = os.listdir(slip_viewer_dir)
        plugin_names = [
            f.split(".py")[0] for f in contents if f.endswith("_frame_plugin.py")
        ]
        self.plugins = {}
        for name in plugin_names:
            self.plugins[name] = SourceFileLoader(
                name, os.path.join(slip_viewer_dir, name + ".py")
            ).load_module()
        if len(plugin_names) > 0:
            print("Loaded plugins: " + ", ".join(plugin_names))

        wx.Frame.__init__(self, *args, **kwds)
        self.settings = rv_settings()

        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(self.sizer)

        # initialization is done in stages as windows are created
        self.pyslip = None
        self.viewer = wx.Panel(self, wx.ID_ANY)
        self.viewer.SetMinSize((640, 640))
        self.viewer.SetBackgroundColour(wx.BLACK)
        self.viewer.ClearBackground()
        self.sizer.Add(self.viewer, 1, wx.EXPAND)

        self.statusbar = self.CreateStatusBar()
        self.settings_frame = None
        self._calibration_frame = None
        self._ring_frame = None
        self._uc_frame = None
        self._score_frame = None
        self._line_frame = None
        self._plugins_frame = dict.fromkeys(self.plugins)
        self.zoom_frame = None
        self.plot_frame = None

        self.metrology_matrices = None

        # Currently displayed image.  XXX Can this be zapped?
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
        self.SetSize((720, 720))

        self.Bind(EVT_EXTERNAL_UPDATE, self.OnExternalUpdate)

        self.Bind(wx.EVT_UPDATE_UI, self.OnUpdateUICalibration, id=self._id_calibration)
        self.Bind(wx.EVT_UPDATE_UI, self.OnUpdateUINext, id=wx.ID_FORWARD)
        self.Bind(wx.EVT_UPDATE_UI, self.OnUpdateUIPrevious, id=wx.ID_BACKWARD)
        self.Bind(wx.EVT_UPDATE_UI, self.OnUpdateUIRing, id=self._id_ring)
        self.Bind(wx.EVT_UPDATE_UI, self.OnUpdateUIUC, id=self._id_uc)
        self.Bind(wx.EVT_UPDATE_UI, self.OnUpdateUIScore, id=self._id_score)
        self.Bind(wx.EVT_UPDATE_UI, self.OnUpdateUILine, id=self._id_line)
        for p in self.plugins:
            self.Bind(
                wx.EVT_UPDATE_UI,
                self.OnUpdateUIPluginWrapper(p),
                id=self._id_plugins[p],
            )

    # consolidate initialization of PySlip object into a single function
    def init_pyslip(self):
        self.set_pyslip(self.viewer)
        self.init_pyslip_presizer()

    def Show(self):
        # Due to the asynchronous nature of X11 on Linux, just showing a frame
        # does not guarantee window creation. The frame calls Raise() so that it
        # will be shown. This addresses an error with PySlip requiring the
        # window to exist before instantiation.
        super().Show()
        self.Raise()

    def setup_toolbar(self):
        XFBaseClass.setup_toolbar(self)

        btn = self.toolbar.AddTool(
            toolId=wx.ID_SAVEAS,
            label="Save As...",
            bitmap=bitmaps.fetch_icon_bitmap("actions", "save_all", 32),
            shortHelp="Save As...",
            kind=wx.ITEM_NORMAL,
        )
        self.Bind(wx.EVT_MENU, self.OnSaveAs, btn)

    # using StaticBox creates a horizontal white bar in Linux
    def make_gui(self, parent):
        parent.sizer = wx.BoxSizer(wx.HORIZONTAL)
        parent.SetSizer(parent.sizer)
        parent.sizer.Add(self.pyslip, 1, wx.EXPAND)

    def init_pyslip_presizer(self):
        self.demo_select_dispatch = {}

        # self.tile_directory = None#"/Users/nksauter/rawdata/demo/insulin_1_001.img"

        # build the GUI
        self.make_gui(self.viewer)

        # finally, bind event to handler
        self.pyslip.Bind(pyslip.EVT_PYSLIP_POSITION, self.handle_position_event)

    def handle_position_event(self, event):
        """Handle a pySlip POSITION event."""

        posn_str = ""
        if event.position:
            (lon, lat) = event.position
            fast_picture, slow_picture = self.pyslip.tiles.lon_lat_to_picture_fast_slow(
                lon, lat
            )

            posn_str = f"Picture:  fast={fast_picture + 0.5:.3f} / slow={slow_picture + 0.5:.3f} pixels."
            coords = self.pyslip.tiles.get_flex_pixel_coordinates(lon, lat)
            if len(coords) >= 2:
                if len(coords) == 3:
                    readout = int(round(coords[2]))
                else:
                    readout = -1
                # the dials convention is that the center of the pixel is 0.5,0.5, so the extent
                # of the 0,0 pixel is from 0.0,0.0 to 1.0,1.0 in pixel space.  However
                # get_flex_pixel_coordinates defines the center of the pixel as 0.0,0.0, so we
                # are a half pixel off
                coords[0] += 0.5
                coords[1] += 0.5
                coords_str = f"fast={coords[1]:.3f} / slow={coords[0]:.3f} pixels"
                if len(coords) == 2:
                    posn_str += " Readout: " + coords_str + "."
                elif readout >= 0:
                    posn_str += " Readout %d: %s." % (readout, coords_str)
                possible_intensity = None
                fi = self.pyslip.tiles.raw_image
                detector = fi.get_detector()
                ifs = (int(coords[1]), int(coords[0]))  # int fast slow
                isf = (int(coords[0]), int(coords[1]))  # int slow fast
                image_data = fi.get_image_data()
                if not isinstance(image_data, tuple):
                    image_data = (image_data,)
                if readout >= 0:
                    if detector[readout].is_coord_valid(ifs):
                        possible_intensity = image_data[readout][isf]

                if possible_intensity is not None:
                    if possible_intensity == 0:
                        format_str = " I=%6.4f"
                        posn_str += format_str % possible_intensity
                    elif possible_intensity == MASK_VAL:
                        posn_str += " I=mask"
                    else:
                        yaya = int(math.ceil(math.log10(abs(possible_intensity))))
                        format_str = " I=%%6.%df" % (max(0, 5 - yaya))
                        posn_str += format_str % possible_intensity

                if (
                    len(coords) > 2 and readout >= 0
                ):  # indicates it's a tiled image in a valid region
                    reso = self.pyslip.tiles.get_resolution(
                        coords[1], coords[0], readout
                    )
                else:
                    reso = self.pyslip.tiles.get_resolution(coords[1], coords[0])

                if reso is not None:
                    posn_str += f" Resolution: {reso:.3f}"

            self.statusbar.SetStatusText(posn_str)
        else:
            self.statusbar.SetStatusText(
                "Click and drag to pan, mouse wheel or double click to zoom"
            )
            # print "event with no position",event
        return

    def init_pyslip_postsizer(self):
        self.pyslip.ZoomToLevel(-2)  # tiles.zoom_level
        self.pyslip.GotoPosition(
            self.pyslip.tiles.get_initial_instrument_centering_within_picture_as_lon_lat()
        )

    def setup_menus(self):
        file_menu = wx.Menu()
        self.mb.Append(file_menu, "File")
        item = file_menu.Append(-1, "Open integration results...")
        self.Bind(wx.EVT_MENU, self.OnLoadIntegration, item)
        item = file_menu.Append(-1, "Open image...")
        self.Bind(wx.EVT_MENU, self.OnLoadFile, item)
        self._actions_menu = wx.Menu()
        self.mb.Append(self._actions_menu, "Actions")
        # item = self._actions_menu.Append(-1, "Change beam center...")
        # self.Bind(wx.EVT_MENU, self.OnChangeBeamCenter, item)
        # item = self._actions_menu.Append(-1, "Reset beam center to header value")
        # self.Bind(wx.EVT_MENU, lambda evt: self.viewer.ResetBeamCenter(), item)
        item = self._actions_menu.Append(-1, "Save As...")
        self.Bind(wx.EVT_MENU, self.OnSaveAs, item)

        # Known wxWidgets/wxPython issue
        # (http://trac.wxwidgets.org/ticket/12394): stock item ID is
        # expected for zero-length text.  Work around by making text
        # contain single space. XXX Placement
        self._id_calibration = wx.NewId()
        item = self._actions_menu.Append(self._id_calibration, " ")
        self.Bind(wx.EVT_MENU, self.OnCalibration, source=item)

        # XXX Placement
        self._id_ring = wx.NewId()
        item = self._actions_menu.Append(self._id_ring, " ")
        self.Bind(wx.EVT_MENU, self.OnRing, source=item)

        # XXX Placement
        self._id_uc = wx.NewId()
        item = self._actions_menu.Append(self._id_uc, " ")
        self.Bind(wx.EVT_MENU, self.OnUC, source=item)

        # XXX Placement
        self._id_score = wx.NewId()
        item = self._actions_menu.Append(self._id_score, " ")
        self.Bind(wx.EVT_MENU, self.OnScore, source=item)

        # XXX Placement
        self._id_line = wx.NewId()
        item = self._actions_menu.Append(self._id_line, " ")
        self.Bind(wx.EVT_MENU, self.OnLine, source=item)

        self._id_plugins = {}
        for p in self.plugins:
            self._id_plugins[p] = wx.NewId()
            item = self._actions_menu.Append(self._id_plugins[p], " ")
            self.Bind(wx.EVT_MENU, self.OnPluginWrapper(p), source=item)

    def has_four_quadrants(self):
        d = self.pyslip.tiles.raw_image.get_detector()
        return len(d) > 1 and len(d.hierarchy()) == 4

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
            if key == str(self.image_chooser.GetClientData(i)):
                return i
        if self.image_chooser.GetCount() >= self.CHOOSER_SIZE:
            self.image_chooser.Delete(0)
        i = self.image_chooser.GetCount()
        if isinstance(file_name_or_data, dict):
            self.image_chooser.Insert(key, i, None)
        elif isinstance(file_name_or_data, chooser_wrapper):
            self.image_chooser.Insert(key, i, file_name_or_data)
        else:
            self.image_chooser.Insert(os.path.basename(key), i, key)
        return i

    def get_beam_center_px(self):
        """
        Get the beam center in pixel coordinates relative to the tile closest to it.
        @return panel_id, beam_center_fast, beam_center_slow. panel_id is the panel the
        returned coordinates are relative to.
        """
        image = self.image_chooser.GetClientData(
            self.image_chooser.GetSelection()
        ).image_set
        detector = image.get_detector()
        beam = image.get_beam()
        if abs(detector[0].get_distance()) == 0:
            return 0.0, 0.0

        # FIXME assumes all detector elements use the same millimeter-to-pixel convention
        try:
            # determine if the beam intersects one of the panels
            panel_id, (x_mm, y_mm) = detector.get_ray_intersection(beam.get_unit_s0())
        except RuntimeError as e:
            if not ("DXTBX_ASSERT(" in str(e) and ") failure" in str(e)):
                # unknown exception from dxtbx
                raise e
            # find the panel whose center is closest to the beam.
            panel_id = 0
            lowest_res = 0
            for p_id, panel in enumerate(detector):
                w, h = panel.get_image_size()
                res = panel.get_resolution_at_pixel(beam, (w // 2, h // 2))
                if res > lowest_res:
                    panel_id = p_id
                    lowest_res = res
            try:
                x_mm, y_mm = detector[panel_id].get_beam_centre(beam.get_unit_s0())
            except RuntimeError:
                # cope with cases like https://github.com/dials/dials/issues/2478
                x_mm, y_mm = detector[panel_id].get_bidirectional_ray_intersection(
                    beam.get_s0()
                )

        beam_pixel_fast, beam_pixel_slow = detector[panel_id].millimeter_to_pixel(
            (x_mm, y_mm)
        )

        return panel_id, beam_pixel_fast, beam_pixel_slow

    def load_image(self, file_name_or_data, get_image_data=None, show_untrusted=False):
        """The load_image() function displays the image from @p
        file_name_or_data.  The chooser is updated appropriately.
        """

        # Due to a bug in wxPython 3.0.2 for Linux
        # http://trac.wxwidgets.org/ticket/16034
        # the creation of the PySlip object is deferred until it is needed and
        # after other windows are created
        if self.pyslip is None:
            self.init_pyslip()
        # The settings dialog is created after PySlip because it may require
        # values from PySlip
        if self.settings_frame is None:
            self.OnShowSettings(None)
        self.Layout()

        if isinstance(file_name_or_data, chooser_wrapper):
            img = rv_image(file_name_or_data)
        else:
            try:
                img = rv_image(file_name_or_data.get_detectorbase())
            except AttributeError:
                img = rv_image(os.path.abspath(file_name_or_data))

        try:
            title = file_name_or_data.full_path
        except AttributeError:
            title = str(file_name_or_data)
        self.SetTitle(title)

        # Update the selection in the chooser.
        i = self.add_file_name_or_data(file_name_or_data)
        self.image_chooser.SetSelection(i)

        self.pyslip.tiles.show_untrusted = show_untrusted
        self.pyslip.tiles.current_brightness = self.settings.brightness
        self.pyslip.tiles.current_color_scheme = self.settings.color_scheme

        self.pyslip.tiles.set_image(
            file_name_or_data=img,
            metrology_matrices=self.metrology_matrices,
            get_image_data=get_image_data,
            show_saturated=(
                self.settings.display == "image"
                and self.settings.image_type == "corrected"
            ),
        )

        # Initialise position zoom level for first image.  XXX Why do we
        # have to coll ZoomToLevel to refresh subsequent images?
        if self._img is None:
            self.init_pyslip_postsizer()
        else:
            self.pyslip.ZoomToLevel(self.pyslip.tiles.zoom_level)

        self._img = img  # XXX

        self.update_statusbar()  # XXX Not always working?
        # self.Layout()

        detector = self.get_detector()
        if abs(detector[0].get_distance()) > 0:

            def map_coords(x, y, p):
                y, x = self.pyslip.tiles.flex_image.tile_readout_to_picture(
                    p, y - 0.5, x - 0.5
                )
                return self.pyslip.tiles.picture_fast_slow_to_map_relative(x, y)

            panel_id, beam_pixel_fast, beam_pixel_slow = self.get_beam_center_px()
            self.beam_center_cross_data = [
                (
                    (
                        map_coords(beam_pixel_fast + 3.0, beam_pixel_slow, panel_id),
                        map_coords(beam_pixel_fast - 3.0, beam_pixel_slow, panel_id),
                    ),
                    {"width": 2, "color": "#0000FFA0", "closed": False},
                ),
                (
                    (
                        map_coords(beam_pixel_fast, beam_pixel_slow + 3.0, panel_id),
                        map_coords(beam_pixel_fast, beam_pixel_slow - 3.0, panel_id),
                    ),
                    {"width": 2, "color": "#0000FFA0", "closed": False},
                ),
            ]

        # Unconditionally delete extra layers--update_settings() will add
        # them back if appropriate.  This also creates the self.*_layer
        # variables.
        if hasattr(self, "beam_layer") and self.beam_layer is not None:
            self.pyslip.DeleteLayer(self.beam_layer, update=False)
        self.beam_layer = None

        if hasattr(self, "spotfinder_layer") and self.spotfinder_layer is not None:
            self.pyslip.DeleteLayer(self.spotfinder_layer)
        self.spotfinder_layer = None

        if hasattr(self, "tile_layer") and self.tile_layer is not None:
            self.pyslip.DeleteLayer(self.tile_layer)
        self.tile_layer = None

        if hasattr(self, "tile_text_layer") and self.tile_text_layer is not None:
            self.pyslip.DeleteLayer(self.tile_text_layer)
        self.tile_text_layer = None

        # if hasattr(self, 'plugins_layer') and hasattr(self.plugins_layer, "__iter__"):
        #   for key in self.plugins_layer:
        #     if self.plugins_layer[key] is not None:
        #       self.pyslip.DeleteLayer(self.plugins_layer[key])
        # self.plugins_layer = {key:None for key in self.plugins}

        self.update_settings()

        # Destroy the calibration frame if it present but unsupported for
        # this image.  XXX Need to do something about the ring tool too
        # when switching between different kinds of images.  XXX Centering
        # is broken when switching between different kinds of images.
        if self._calibration_frame and not self.has_four_quadrants():
            self.OnCalibration(None)

    def get_detector(self):
        return self.pyslip.tiles.raw_image.get_detector()

    def get_beam(self):
        return self.pyslip.tiles.raw_image.get_beam()

    def get_key(self, file_name_or_data):
        """This overridden get_key() function returns the key of @p file_name_or_data
        if it's an DetectorImageBase object.  Otherwise it returns the super class's
        key
        """
        from iotbx.detectors.detectorbase import DetectorImageBase

        if isinstance(file_name_or_data, DetectorImageBase):
            return file_name_or_data.filename
        elif isinstance(file_name_or_data, chooser_wrapper):
            return str(file_name_or_data)
        else:
            return super().get_key(file_name_or_data)

    def update_settings(self, layout=True):
        raise NotImplementedError()

    def OnCalibration(self, event):
        if not self._calibration_frame:
            self._calibration_frame = SBSettingsFrame(
                self, wx.ID_ANY, "Quadrant calibration", style=wx.CAPTION | wx.CLOSE_BOX
            )
            self._calibration_frame.Show()
            self._calibration_frame.Raise()
        else:
            self._calibration_frame.Destroy()

    def OnRing(self, event):
        from .ring_frame import RingSettingsFrame

        if not self._ring_frame:
            self._ring_frame = RingSettingsFrame(
                self, wx.ID_ANY, "Ring tool", style=wx.CAPTION | wx.CLOSE_BOX
            )
            self._ring_frame.Show()
            self._ring_frame.Raise()
        else:
            self._ring_frame.Destroy()

    def OnUC(self, event):
        from .uc_frame import UCSettingsFrame

        if not self._uc_frame:
            self._uc_frame = UCSettingsFrame(
                self, wx.ID_ANY, "Unit cell tool", style=wx.CAPTION | wx.CLOSE_BOX
            )
            self._uc_frame.Show()
            self._uc_frame.Raise()
        else:
            self._uc_frame.Destroy()

    def OnScore(self, event):
        from .score_frame import ScoreSettingsFrame

        if not self._score_frame:
            self._score_frame = ScoreSettingsFrame(
                self, wx.ID_ANY, "Score tool", style=wx.CAPTION | wx.CLOSE_BOX
            )
            self._score_frame.Show()
            self._score_frame.Raise()
        else:
            self._score_frame.Destroy()

    def OnLine(self, event):
        from .line_frame import LineSettingsFrame

        if not self._line_frame:
            self._line_frame = LineSettingsFrame(
                self,
                wx.ID_ANY,
                "Line tool",
                style=wx.CAPTION | wx.CLOSE_BOX | wx.RESIZE_BORDER,
            )
            self._line_frame.Show()
            self._line_frame.Raise()
        else:
            self._line_frame.Destroy()

    def OnPluginWrapper(self, p):
        def OnPlugin(event):
            if not self._plugins_frame[p]:
                helper = self.plugins[p].PluginHelper
                self._plugins_frame[p] = helper._plugin_settings_frame(
                    self,
                    wx.ID_ANY,
                    helper._plugin_title,
                    style=wx.CAPTION | wx.CLOSE_BOX,
                )
                self._plugins_frame[p].Show()
                self._plugins_frame[p].Raise()
            else:
                self._plugins_frame[p].Destroy()

        return OnPlugin

    def OnUpdateUICalibration(self, event):
        # If quadrant calibration is not supported for this image, disable
        # the corresponding menu item.  Toggle the menu item text
        # depending on the state of the tool.

        if self.has_four_quadrants():
            event.Enable(True)
            if self._calibration_frame:
                event.SetText("Hide quadrant calibration")
            else:
                event.SetText("Show quadrant calibration")
        else:
            event.Enable(False)
            event.SetText("Show quadrant calibration")

    def OnUpdateUINext(self, event):
        # Enable/disable the "Next" button based on the image's position
        # in the list.

        event.Enable(
            self.image_chooser.GetSelection() + 1 < self.image_chooser.GetCount()
        )

    def OnUpdateUIPrevious(self, event):
        # Enable/disable the "Previous" button based on the image's
        # position in the list.

        event.Enable(self.image_chooser.GetSelection() >= 1)

    def OnUpdateUIRing(self, event):
        # Toggle the menu item text depending on the state of the tool.

        if self._ring_frame:
            event.SetText("Hide ring tool")
        else:
            event.SetText("Show ring tool")

    def OnUpdateUIUC(self, event):
        # Toggle the menu item text depending on the state of the tool.

        if self._uc_frame:
            event.SetText("Hide unit cell tool")
        else:
            event.SetText("Show unit cell tool")

    def OnUpdateUIScore(self, event):
        # Toggle the menu item text depending on the state of the tool.

        if self._score_frame:
            event.SetText("Hide score tool")
        else:
            event.SetText("Show score tool")

    def OnUpdateUILine(self, event):
        # Toggle the menu item text depending on the state of the tool.

        if self._line_frame:
            event.SetText("Hide line tool")
        else:
            event.SetText("Show line tool")

    def OnUpdateUIPluginWrapper(self, p):
        def OnUpdateUIPlugin(event):
            # Toggle the menu item text depending on the state of the tool.

            helper = self.plugins[p].PluginHelper
            if self._plugins_frame[p]:
                event.SetText(helper._plugin_hide_text)
            else:
                event.SetText(helper._plugin_show_text)

        return OnUpdateUIPlugin

    def OnSaveAs(self, event):
        ### XXX TODO: Save overlays
        ### XXX TODO: Fix bug where multi-asic images are slightly cropped due to transformation error'

        import PIL.Image as Image

        dialog = wx.FileDialog(
            self,
            defaultDir="",
            message="Save PNG or PDF file",
            style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT,
            wildcard="PNG file (*.png)|*.png|PDF file (*.pdf)|*.pdf",
        )
        if dialog.ShowModal() != wx.ID_OK:
            return

        file_name = dialog.GetPath()
        if file_name == "":
            return

        self.update_statusbar("Writing " + file_name + "...")
        if dialog.GetFilterIndex() == 0:
            # XXX Copied from tile_generation.py; all its disclaimers
            # apply.
            raw_img = self.pyslip.tiles.raw_image
            detector = raw_img.get_detector()
            data = raw_img.get_image_data()
            if not isinstance(data, tuple):  # XXX should not need this test
                data = (data,)
            from .tile_generation import get_flex_image_multipanel

            flex_img = get_flex_image_multipanel(
                brightness=self.settings.brightness / 100,
                detector=detector,
                image_data=data,
                beam=raw_img.get_beam(),
            )

            if flex_img.supports_rotated_tiles_antialiasing_recommended:
                currentZoom = self.pyslip.level
                self.pyslip.tiles.UseLevel(0)  # 1:1 zoom level

                try:
                    x, y, width, height = self._img._raw.bounding_box_mm()
                    x1, y1 = self._img._raw.detector_coords_as_image_coords(x, y)
                    x2, y2 = self._img._raw.detector_coords_as_image_coords(
                        x + width, y + height
                    )
                except AttributeError:
                    x1 = min(
                        [
                            p.get_pixel_lab_coord(c)[0] / p.get_pixel_size()[0]
                            for p in detector
                            for c in [
                                (0, 0),
                                (0, p.get_image_size()[1]),
                                (p.get_image_size()[0], 0),
                                (p.get_image_size()[0], p.get_image_size()[1]),
                            ]
                        ]
                    )
                    y1 = min(
                        [
                            p.get_pixel_lab_coord(c)[1] / p.get_pixel_size()[1]
                            for p in detector
                            for c in [
                                (0, 0),
                                (0, p.get_image_size()[1]),
                                (p.get_image_size()[0], 0),
                                (p.get_image_size()[0], p.get_image_size()[1]),
                            ]
                        ]
                    )
                    x2 = max(
                        [
                            p.get_pixel_lab_coord(c)[0] / p.get_pixel_size()[0]
                            for p in detector
                            for c in [
                                (0, 0),
                                (0, p.get_image_size()[1]),
                                (p.get_image_size()[0], 0),
                                (p.get_image_size()[0], p.get_image_size()[1]),
                            ]
                        ]
                    )
                    y2 = max(
                        [
                            p.get_pixel_lab_coord(c)[1] / p.get_pixel_size()[1]
                            for p in detector
                            for c in [
                                (0, 0),
                                (0, p.get_image_size()[1]),
                                (p.get_image_size()[0], 0),
                                (p.get_image_size()[0], p.get_image_size()[1]),
                            ]
                        ]
                    )

                # Map > View - determine layout in X direction
                x_offset = x1
                start_x_tile = int(math.floor(x_offset / self.pyslip.tile_size_x))
                stop_x_tile = (
                    x2 + self.pyslip.tile_size_x - 1
                ) / self.pyslip.tile_size_x
                stop_x_tile = int(stop_x_tile)
                col_list = range(start_x_tile, stop_x_tile)
                x_pix = start_x_tile * self.pyslip.tile_size_y - x_offset

                y_offset = y1
                start_y_tile = int(math.floor(y_offset / self.pyslip.tile_size_y))
                stop_y_tile = (
                    y2 + self.pyslip.tile_size_y - 1
                ) / self.pyslip.tile_size_y
                stop_y_tile = int(stop_y_tile)
                row_list = range(start_y_tile, stop_y_tile)
                y_pix_start = start_y_tile * self.pyslip.tile_size_y - y_offset

                bitmap = wx.Bitmap(int(x2 - x1), int(y2 - y1))
                dc = wx.MemoryDC()
                dc.SelectObject(bitmap)

                # start pasting tiles
                for x in col_list:
                    y_pix = y_pix_start
                    for y in row_list:
                        dc.DrawBitmap(
                            self.pyslip.tiles.GetTile(x, y),
                            int(x_pix),
                            int(y_pix),
                            False,
                        )
                        y_pix += self.pyslip.tile_size_y
                    x_pix += self.pyslip.tile_size_x

                dc.SelectObject(wx.NullBitmap)

                wximg = bitmap.ConvertToImage()
                imageout = Image.new("RGB", (wximg.GetWidth(), wximg.GetHeight()))
                imageout.frombytes(bytes(wximg.GetData()))

                self.pyslip.tiles.UseLevel(currentZoom)

            else:  # write the image out at full resolution
                flex_img.setWindow(0.0, 0.0, 1)
                flex_img.spot_convention(0)
                flex_img.adjust(color_scheme=self.settings.color_scheme)
                flex_img.prep_string()
                data_string = flex_img.as_bytes()
                imageout = Image.frombytes(
                    "RGB", (flex_img.ex_size2(), flex_img.ex_size1()), data_string
                )

            with open(file_name, "wb") as fh:
                imageout.save(fh, "PNG")

        elif dialog.GetFilterIndex() == 1:
            from reportlab.lib.units import inch
            from reportlab.pdfgen import canvas

            # Dots per inch in PDF output, and fudge factor to not make
            # fine features impossibly small.  XXX The fudge factor should
            # go.
            DPI = 72
            LINE_WIDTH_FACTOR = 0.6

            # XXX Copied from tile_generation.py; all its disclaimers
            # apply.
            raw_img = self.pyslip.tiles.raw_image
            detector = raw_img.get_detector()
            data = raw_img.get_image_data()
            if not isinstance(data, tuple):  # XXX should not need this test
                data = (data,)
            from .tile_generation import get_flex_image_multipanel

            flex_img = get_flex_image_multipanel(
                brightness=self.settings.brightness / 100,
                detector=detector,
                image_data=data,
                beam=raw_img.get_beam(),
            )

            flex_img.setWindow(0, 0, 1)
            flex_img.adjust(color_scheme=self.settings.color_scheme)
            flex_img.prep_string()

            # XXX Order of size1/size2 correct?
            pdf_size = (flex_img.size2() * inch / DPI, flex_img.size1() * inch / DPI)
            pdf_canvas = canvas.Canvas(filename=file_name, pagesize=pdf_size)
            pil_img = Image.frombytes(
                "RGB", (flex_img.size2(), flex_img.size1()), flex_img.as_bytes()
            )

            pdf_canvas.drawInlineImage(
                pil_img, 0, 0, width=pdf_size[0], height=pdf_size[1]
            )

            for layer_id in self.pyslip.layer_z_order:
                layer = self.pyslip.layer_mapping[layer_id]

                # XXX This would probably be more elegant if these were
                # functions in some layer class.  Note repeated roundabout
                # way (via a wx.Pen object) to extract RGB values from the
                # colour parameter.
                if layer.type == self.pyslip.TypeEllipse:
                    for (
                        p,
                        place,
                        width,
                        colour,
                        closed,
                        filled,
                        fillcolour,
                        x_off,
                        y_off,
                        pdata,
                    ) in layer.data:
                        if layer.map_rel:
                            pp = []
                            for pelement in p:
                                fs = (
                                    self.pyslip.tiles.map_relative_to_picture_fast_slow(
                                        pelement[0], pelement[1]
                                    )
                                )
                                pp.append(
                                    (
                                        fs[0] * inch / DPI,
                                        pdf_size[1] - fs[1] * inch / DPI,
                                    )
                                )
                            ellipse_center = (pp[0][0], pp[0][1])
                            major = (
                                pp[1][0] - ellipse_center[0],
                                pp[1][1] - ellipse_center[1],
                            )
                            minor = (
                                pp[2][0] - ellipse_center[0],
                                pp[2][1] - ellipse_center[1],
                            )
                        else:
                            raise NotImplementedError(
                                "PDF output in view-relative coordinates not implemented"
                            )

                        pen = wx.Pen(colour)
                        pdf_canvas.setLineWidth(width * LINE_WIDTH_FACTOR)
                        pdf_canvas.setStrokeColorRGB(
                            pen.Colour.Red() / 255,
                            pen.Colour.Green() / 255,
                            pen.Colour.Blue() / 255,
                        )

                        angle = math.atan2(major.elems[1], major.elems[0])
                        r_major = math.hypot(major.elems[0], major.elems[1])
                        r_minor = math.hypot(minor.elems[0], minor.elems[1])

                        pdf_canvas.saveState()
                        pdf_canvas.translate(
                            ellipse_center.elems[0], ellipse_center.elems[1]
                        )

                        pdf_canvas.rotate(math.degrees(angle))
                        pdf_canvas.ellipse(-r_major, -r_minor, r_major, r_minor)
                        pdf_canvas.restoreState()

                elif layer.type == self.pyslip.TypeImage:
                    raise NotImplementedError(
                        "PDF output of image layers not implemented"
                    )

                elif layer.type == self.pyslip.TypePoint:
                    for (
                        lon,
                        lat,
                        place,
                        radius,
                        colour,
                        x_off,
                        y_off,
                        pdata,
                    ) in layer.data:
                        if layer.map_rel:
                            fs = self.pyslip.tiles.map_relative_to_picture_fast_slow(
                                lon, lat
                            )
                        else:
                            raise NotImplementedError(
                                "PDF output in view-relative coordinates not implemented"
                            )

                        pt = (fs[0] * inch / DPI, pdf_size[1] - fs[1] * inch / DPI)
                        pen = wx.Pen(colour)
                        pdf_canvas.setLineWidth(radius)
                        pdf_canvas.setStrokeColorRGB(
                            pen.Colour.Red() / 255,
                            pen.Colour.Green() / 255,
                            pen.Colour.Blue() / 255,
                        )
                        pdf_canvas.circle(pt[0], pt[1], 0.5 * radius * inch / DPI)

                elif layer.type == self.pyslip.TypePolygon:
                    for (
                        p,
                        place,
                        width,
                        colour,
                        closed,
                        filled,
                        fillcolour,
                        x_off,
                        y_off,
                        pdata,
                    ) in layer.data:
                        path = pdf_canvas.beginPath()
                        for i, pp in enumerate(p):
                            if layer.map_rel:
                                fs = (
                                    self.pyslip.tiles.map_relative_to_picture_fast_slow(
                                        pp[0], pp[1]
                                    )
                                )
                            else:
                                raise NotImplementedError(
                                    "PDF output in view-relative coordinates not implemented"
                                )

                            pt = (fs[0] * inch / DPI, pdf_size[1] - fs[1] * inch / DPI)
                            if i == 0:
                                path.moveTo(pt[0], pt[1])
                            else:
                                path.lineTo(pt[0], pt[1])

                        if closed:
                            path.close()

                        pen = wx.Pen(colour)
                        pdf_canvas.setFillColorRGB(
                            pen.Colour.Red() / 255,
                            pen.Colour.Green() / 255,
                            pen.Colour.Blue() / 255,
                        )
                        pdf_canvas.setLineWidth(width * LINE_WIDTH_FACTOR)
                        pdf_canvas.setStrokeColorRGB(
                            pen.Colour.Red() / 255,
                            pen.Colour.Green() / 255,
                            pen.Colour.Blue() / 255,
                        )
                        pdf_canvas.drawPath(path, fill=filled)

                elif layer.type == self.pyslip.TypeText:
                    for (
                        lon,
                        lat,
                        tdata,
                        placement,
                        radius,
                        colour,
                        textcolour,
                        fontname,
                        fontsize,
                        offset_x,
                        offset_y,
                        data,
                    ) in layer.data:
                        if placement != "cc":
                            print(
                                Warning(
                                    "Only centered placement available when drawing text on pdf"
                                )
                            )
                        if layer.map_rel:
                            fs = self.pyslip.tiles.map_relative_to_picture_fast_slow(
                                lon, lat
                            )
                        else:
                            raise NotImplementedError(
                                "PDF output in view-relative coordinates not implemented"
                            )
                        from reportlab.pdfbase.pdfmetrics import stringWidth

                        scale = 5  # XXX this scaleup by 5 is arbitrary!
                        try:
                            w = stringWidth(tdata, fontname, fontsize * scale)
                        except KeyError:
                            fontname = "Helvetica"
                            w = stringWidth(tdata, fontname, fontsize * scale)
                        if fs[0] - (w / 2) < 0:  # handle text falling off the left side
                            txt = pdf_canvas.beginText(x=0, y=fs[1])
                        else:
                            txt = pdf_canvas.beginText(x=fs[0] - (w / 2), y=fs[1])
                        txt.setFont(fontname, fontsize * scale)
                        if isinstance(textcolour, wx.Colour):
                            textcolour = tuple(textcolour)
                        txt.setFillColor(textcolour)
                        txt.setStrokeColor(textcolour)
                        txt.textLine(tdata)
                        pdf_canvas.drawText(txt)

            pdf_canvas.save()

        self.update_statusbar("Writing " + file_name + "..." + " Done.")
