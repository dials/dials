# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1

"""
Visualise the reference profiles from integration.

Examples::

  dials.integrate refined.expt refined.refl debug.reference.output=True
  dials.reference_profile_viewer reference_profiles.pickle
"""

from __future__ import annotations

import copy
import os
import pickle

import matplotlib
import numpy
import wx

matplotlib.use("WXAgg")

from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas
from matplotlib.backends.backend_wxagg import (
    NavigationToolbar2WxAgg as NavigationToolbar,
)
from matplotlib.figure import Figure

import dials.util
import dials.util.log
from dials.util.options import ArgumentParser


class ProfilesFrame(wx.Frame):
    """The main frame of the application"""

    title = "Reference profile viewer"

    def __init__(self, profiles):
        wx.Frame.__init__(self, None, -1, self.title)

        self.profiles = profiles
        self.data = self.profiles.get_profiles(experiment=0, block=0)

        self.create_menu()
        self.create_status_bar()
        self.create_main_panel()
        self.draw_figure()

    def create_menu(self):
        self.menubar = wx.MenuBar()

        menu_file = wx.Menu()
        m_expt = menu_file.Append(-1, "&Save plot\tCtrl-S", "Save plot to file")
        self.Bind(wx.EVT_MENU, self.on_save_plot, m_expt)
        menu_file.AppendSeparator()
        m_exit = menu_file.Append(-1, "E&xit\tCtrl-X", "Exit")
        self.Bind(wx.EVT_MENU, self.on_exit, m_exit)

        menu_help = wx.Menu()
        m_about = menu_help.Append(-1, "&About\tF1", "About the program")
        self.Bind(wx.EVT_MENU, self.on_about, m_about)

        self.menubar.Append(menu_file, "&File")
        self.menubar.Append(menu_help, "&Help")
        self.SetMenuBar(self.menubar)

    def create_main_panel(self):
        """Creates the main panel with all the controls on it:
        * mpl canvas
        * mpl navigation toolbar
        * Control panel for interaction
        """
        self.panel = wx.Panel(self)

        # Create the mpl Figure and FigCanvas objects.
        # 7x5 inches, 100 dots-per-inch
        self.dpi = 100
        self.fig = Figure((7.0, 5.0), dpi=self.dpi)
        self.canvas = FigCanvas(self.panel, -1, self.fig)

        self.set_axes()

        self.expt_selection_label = wx.StaticText(self.panel, -1, "Experiment ID: ")
        self.expt_selection = wx.SpinCtrl(
            self.panel, -1, "0", max=self.profiles.get_n_experiments() - 1
        )
        self.Bind(wx.EVT_SPINCTRL, self.on_spin_expt, self.expt_selection)

        self.block_selection_label = wx.StaticText(self.panel, -1, "Block: ")
        self.block_selection = wx.SpinCtrl(
            self.panel, -1, "0", max=self.profiles.get_n_blocks(experiment=0) - 1
        )
        self.Bind(wx.EVT_SPINCTRL, self.on_spin_block, self.block_selection)

        self.mask_checkbox = wx.CheckBox(self.panel, -1, "Mask")
        self.mask_checkbox.SetValue(True)
        self.Bind(wx.EVT_CHECKBOX, self.redraw_on_event, self.mask_checkbox)

        self.cmap_choice_label = wx.StaticText(self.panel, -1, "Colourmap: ")
        self.cmap_choice = wx.Choice(
            self.panel, -1, choices=["viridis", "viridis_r", "gray", "gray_r"]
        )
        self.cmap_choice.SetSelection(0)
        self.Bind(wx.EVT_CHOICE, self.redraw_on_event, self.cmap_choice)

        # Create the navigation toolbar, tied to the canvas
        self.toolbar = NavigationToolbar(self.canvas)

        # Layout with box sizers
        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.vbox.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
        self.vbox.Add(self.toolbar, 0, wx.EXPAND)
        self.vbox.AddSpacer(10)

        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        flags = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL
        self.hbox.Add(self.expt_selection_label, 0, flag=flags)
        self.hbox.Add(self.expt_selection, 0, border=3, flag=flags)
        self.hbox.Add(self.block_selection_label, 0, flag=flags)
        self.hbox.Add(self.block_selection, 0, border=3, flag=flags)
        self.hbox.Add(self.mask_checkbox, 0, border=3, flag=flags)
        self.hbox.Add(self.cmap_choice_label, 0, flag=flags)
        self.hbox.Add(self.cmap_choice, 0, border=3, flag=flags)

        self.vbox.Add(self.hbox, 0, flag=wx.ALIGN_LEFT | wx.TOP)

        self.panel.SetSizer(self.vbox)
        self.vbox.Fit(self)

    def set_axes(self):
        subplots = [e["subplot"] for e in self.data]
        r, c = list(zip(*subplots))
        nrows = max(r) + 1
        ncols = max(c) + 1
        self.axes = self.fig.subplots(nrows, ncols, sharex=True, sharey=True)

    def create_status_bar(self):
        self.statusbar = self.CreateStatusBar()

    def draw_figure(self):
        """Redraws the figure"""

        final_row_index = self.axes.shape[0] - 1
        for profile in self.data:
            subplot = profile["subplot"]
            ax = self.axes[subplot]
            ax.clear()

            # For now, let's just sum down the first axis. Profiles are stored
            # in ez, ey, ex order, where ex is orthogonal to s1 and s0, ey is
            # orthogonal to s1 and ex, and ez is the axis that is dependent on
            # the direction through the Ewald sphere
            vals2D = profile["data"].sum(axis=0)
            cmap = copy.copy(
                matplotlib.cm.get_cmap(
                    self.cmap_choice.GetString(self.cmap_choice.GetSelection())
                )
            )

            # If any X, Y position is masked down the summed Z stack then mask
            # it in the final image.
            if self.mask_checkbox.IsChecked():
                mask2D = (profile["mask"] - 1).sum(axis=0)
                mask2D = mask2D != 0
                vals2D[mask2D] = numpy.nan
            cmap.set_bad(color="red")
            ax.imshow(vals2D, cmap=cmap)

            if subplot[0] == final_row_index:
                ax.set_xlabel(f"X (px): {profile['coord'][0]:.1f}")
            if subplot[1] == 0:
                ax.set_ylabel(f"Y (px): {profile['coord'][1]:.1f}")

        self.fig.suptitle(f"Block Z (im): {profile['coord'][2]:.1f}")
        self.canvas.draw()

    def on_spin_expt(self, event):
        self.expt_selection.Disable()
        exp_id = self.expt_selection.GetValue()
        self.data = self.profiles.get_profiles(experiment=exp_id, block=0)
        self.block_selection.SetValue(0)
        self.block_selection.SetMax(self.profiles.get_n_blocks(exp_id) - 1)
        self.draw_figure()
        self.expt_selection.Enable()

    def on_spin_block(self, event):
        self.block_selection.Disable()
        self.data = self.profiles.get_profiles(
            experiment=self.expt_selection.GetValue(),
            block=self.block_selection.GetValue(),
        )
        self.draw_figure()
        self.block_selection.Enable()

    def redraw_on_event(self, event):
        self.draw_figure()

    def on_save_plot(self, event):
        file_choices = "PNG (*.png)|*.png"

        dlg = wx.FileDialog(
            self,
            message="Save plot as...",
            defaultDir=os.getcwd(),
            defaultFile="reference_profiles.png",
            wildcard=file_choices,
            style=wx.SAVE,
        )

        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            self.canvas.print_figure(path, dpi=self.dpi)
            self.flash_status_message("Saved to %s" % path)

    def on_exit(self, event):
        self.Destroy()

    def on_about(self, event):
        msg = """A reference profile viewer for DIALS:

         * The reference profiles for the selected block of
           images (Z) are displayed along with their X, Y
           pixel positions
         * The display sums slices of the Kabsch space shoebox
           down the e3 direction, hence the view is the
           integration of the profile as it passes through
           the Ewald sphere
         * If any of the slices contains a masked pixel, the
           equivalent pixel in the summed image is also
           showed as masked
         * Save the plot to a file using the File menu
        """
        dlg = wx.MessageDialog(self, msg, "About", wx.OK)
        dlg.ShowModal()
        dlg.Destroy()

    def flash_status_message(self, msg, flash_len_ms=1500):
        self.statusbar.SetStatusText(msg)
        self.timeroff = wx.Timer(self)
        self.Bind(wx.EVT_TIMER, self.on_flash_status_off, self.timeroff)
        self.timeroff.Start(flash_len_ms, oneShot=True)

    def on_flash_status_off(self, event):
        self.statusbar.SetStatusText("")


class ProfileStore:
    def __init__(self, filename):
        with open(filename, "rb") as f:
            data = pickle.load(f)

        self._experiment_data = []
        for d in data:
            self._experiment_data.append(self._process_experiment(d))

    def _process_experiment(self, data):
        # Reorganise data by Z block
        z_blocks = {e["coord"][2] for e in data}
        z_blocks = dict.fromkeys(z_blocks)
        for model in data:
            z_coord = model["coord"][2]
            if z_blocks[z_coord] is None:
                z_blocks[z_coord] = [
                    model,
                ]
            else:
                z_blocks[z_coord].append(model)

        # Process within Z blocks
        for k, v in z_blocks.items():
            z_blocks[k] = self._process_block(v)

        # Convert to list
        keys = sorted(z_blocks.keys())
        return [z_blocks[k] for k in keys]

    @staticmethod
    def _process_block(block):
        x_coords = sorted({e["coord"][0] for e in block})
        y_coords = sorted({e["coord"][1] for e in block})

        for profile in block:
            x, y, _ = profile["coord"]
            i = x_coords.index(x)
            j = y_coords.index(y)
            profile["subplot"] = (j, i)  # row (Y), then column (X)
            profile["data"] = profile["data"].as_numpy_array()
            profile["mask"] = profile["mask"].as_numpy_array()

        return block

    def get_n_experiments(self):
        return len(self._experiment_data)

    def get_n_blocks(self, experiment):
        return len(self._experiment_data[experiment])

    def get_profiles(self, experiment, block):
        return self._experiment_data[experiment][block]


@dials.util.show_mail_on_error()
def run():
    dials.util.log.print_banner()
    usage = "dials.reference_profile_viewer [options] reference_profiles.pickle"
    parser = ArgumentParser(
        usage=usage,
        read_reflections=False,
        read_experiments=False,
        check_format=False,
        epilog=__doc__,
    )

    params, options, args = parser.parse_args(
        show_diff_phil=False, return_unhandled=True
    )
    if len(args) != 1:
        parser.print_help()
        exit(0)
    filename = args[0]
    assert os.path.isfile(filename), filename

    # Load data
    data = ProfileStore(filename)

    # Start viewer
    show_reference_profile_viewer(data, params)


def show_reference_profile_viewer(data, params):
    app = wx.App()
    app.frame = ProfilesFrame(data)
    app.frame.Show()
    app.MainLoop()


if __name__ == "__main__":
    run()
