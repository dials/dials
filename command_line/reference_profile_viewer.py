# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1

"""
This demo demonstrates how to embed a matplotlib (mpl) plot
into a wxPython GUI application, including:

* Using the navigation toolbar
* Adding data to the plot
* Dynamically modifying the plot's properties
* Processing mpl events
* Saving the plot to a file from a menu

The main goal is to serve as a basis for developing rich wx GUI
applications featuring mpl plots (using the mpl OO API).

Eli Bendersky (eliben@gmail.com)
License: this code is in the public domain
"""
import os
import pickle

import matplotlib
import wx

matplotlib.use("WXAgg")

from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas
from matplotlib.backends.backend_wxagg import (
    NavigationToolbar2WxAgg as NavigationToolbar,
)
from matplotlib.figure import Figure

import libtbx.phil

import dials.util
import dials.util.log
from dials.util.options import OptionParser

try:
    from typing import List
except ImportError:
    pass


# Define the master PHIL scope for this program.
phil_scope = libtbx.phil.parse(
    """
    bool_parameter = False
      .type = bool
    integer_parameter = 0
      .type = int
    """
)


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

        self.textbox.SetValue("experiment: 0 block: 0")
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
        #
        self.dpi = 100
        self.fig = Figure((7.0, 5.0), dpi=self.dpi)
        self.canvas = FigCanvas(self.panel, -1, self.fig)

        self.set_axes()

        # Bind the 'pick' event for clicking on one of the values
        #
        self.canvas.mpl_connect("pick_event", self.on_pick)

        self.textbox = wx.TextCtrl(
            self.panel, size=(200, -1), style=wx.TE_PROCESS_ENTER
        )
        self.Bind(wx.EVT_TEXT_ENTER, self.on_text_enter, self.textbox)

        self.drawbutton = wx.Button(self.panel, -1, "Draw!")
        self.Bind(wx.EVT_BUTTON, self.on_draw_button, self.drawbutton)

        self.cb_grid = wx.CheckBox(self.panel, -1, "Show Grid", style=wx.ALIGN_RIGHT)
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_grid, self.cb_grid)

        self.slider_label = wx.StaticText(self.panel, -1, "Bar width (%): ")
        self.slider_width = wx.Slider(
            self.panel,
            -1,
            value=20,
            minValue=1,
            maxValue=100,
            style=wx.SL_AUTOTICKS | wx.SL_LABELS,
        )
        self.slider_width.SetTickFreq(10)
        self.Bind(
            wx.EVT_COMMAND_SCROLL_THUMBTRACK, self.on_slider_width, self.slider_width
        )

        # Create the navigation toolbar, tied to the canvas
        #
        self.toolbar = NavigationToolbar(self.canvas)

        #
        # Layout with box sizers
        #

        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.vbox.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
        self.vbox.Add(self.toolbar, 0, wx.EXPAND)
        self.vbox.AddSpacer(10)

        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        flags = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL
        self.hbox.Add(self.textbox, 0, border=3, flag=flags)
        self.hbox.Add(self.drawbutton, 0, border=3, flag=flags)
        self.hbox.Add(self.cb_grid, 0, border=3, flag=flags)
        self.hbox.AddSpacer(30)
        self.hbox.Add(self.slider_label, 0, flag=flags)
        self.hbox.Add(self.slider_width, 0, border=3, flag=flags)

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
        # str = self.textbox.GetValue()
        # self.data = list(map(int, str.split()))
        # x = range(len(self.data))

        final_row_index = self.axes.shape[0] - 1
        for profile in self.data:
            subplot = profile["subplot"]
            ax = self.axes[subplot]
            ax.clear()

            # for now, let's just sum down the first axis. Profiles are stored
            # in ez, ey, ex order, where ex is orthogonal to s1 and s0, ey is
            # orthogonal to s1 and ex, and ez is the axis that is dependent on
            # the direction through the Ewald sphere
            vals2D = profile["data"].sum(axis=0)
            ax.imshow(vals2D)

            if subplot[0] == final_row_index:
                ax.set_xlabel(f"X: {profile['coord'][0]:.1f}")
            if subplot[1] == 0:
                ax.set_ylabel(f"Y: {profile['coord'][1]:.1f}")

        self.fig.suptitle(f"Block Z: {profile['coord'][2]:.1f}")
        self.fig.tight_layout()

        self.canvas.draw()

    def on_cb_grid(self, event):
        self.draw_figure()

    def on_slider_width(self, event):
        self.draw_figure()

    def on_draw_button(self, event):
        self.draw_figure()

    def on_pick(self, event):
        # The event received here is of the type
        # matplotlib.backend_bases.PickEvent
        #
        # It carries lots of information, of which we're using
        # only a small amount here.
        #
        box_points = event.artist.get_bbox().get_points()
        msg = "You've clicked on a bar with coords:\n %s" % box_points

        dlg = wx.MessageDialog(self, msg, "Click!", wx.OK | wx.ICON_INFORMATION)

        dlg.ShowModal()
        dlg.Destroy()

    def on_text_enter(self, event):
        self.draw_figure()

    def on_save_plot(self, event):
        file_choices = "PNG (*.png)|*.png"

        dlg = wx.FileDialog(
            self,
            message="Save plot as...",
            defaultDir=os.getcwd(),
            defaultFile="plot.png",
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

         * Use the matplotlib navigation bar
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
        z_blocks = set(e["coord"][2] for e in data)
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
        x_coords = sorted(set(e["coord"][0] for e in block))
        y_coords = sorted(set(e["coord"][1] for e in block))

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
def run(args: List[str] = None, phil: libtbx.phil.scope = phil_scope) -> None:
    """
    Check command-line input and call other functions to do the legwork.
    Run the script, parsing arguments found in 'args' and using the PHIL scope
    defined in 'phil'.
    Try to keep this function minimal, defining only what is necessary to run
    the program from the command line.
    Args:
        args: The arguments supplied by the user (default: sys.argv[1:])
        phil: The PHIL scope definition (default: phil_scope, the master PHIL scope
        for this program).
    """
    dials.util.log.print_banner()
    usage = "dials.referemce_profile_viewer [options] reference_profiles.pickle"
    parser = OptionParser(
        usage=usage,
        phil=phil,
        read_reflections=False,
        read_experiments=False,
        check_format=False,
        epilog=__doc__,
    )

    params, options, args = parser.parse_args(
        args=args, show_diff_phil=False, return_unhandled=True
    )
    if len(args) != 1:
        parser.print_help()
        exit(0)
    filename = args[0]
    assert os.path.isfile(filename), filename

    data = ProfileStore(filename)

    # Do whatever this program is supposed to do.
    show_reference_profile_viewer(data, params)


def show_reference_profile_viewer(data, params):
    app = wx.App()
    app.frame = ProfilesFrame(data)
    app.frame.Show()
    app.MainLoop()


# Keep this minimal. Calling run() should do exactly the same thing as running this
if __name__ == "__main__":
    run()
