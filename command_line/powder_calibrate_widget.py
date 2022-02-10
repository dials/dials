# LIBTBX_SET_DISPATCHER_NAME dials.powder_calibrate_widget

"""
Calibrate electron powder geometry using X-ray powder diffraction tools.
    `pyFAI` is a well established X-ray powder diffraction tool.
https://doi.org/10.1107/S1600576715004306

The matplotlib widget requires user eyes to give an initial guess of the beam center.
From the initial guess pyFAI geometry calibration can be used in the usual manner
(ie. just like it used for X-rays).

Usage:
    from dials.command_line.powder_calibrate_widget import DialsParams, PowderCalibrator

    Al_data = DialsParams(expt_file="imported.expt")
    calibrator = PowderCalibrator(Al_data, 'Al')
    calibrator.calibrate_with_calibrant(verbose=True)

Or from command line:
    1. use the widget to generate a starting geometry. This step can be skipped if
    the starting geometry is pretty close.
    > dials.powder_calibrate_widget data=poor_geom.expt standard="Al" calibrated_geom="eyeballed.expt"

    2. start fine calibration from a saved eyeballed geometry
    > dials.powder_calibrate_widget data=eyeballed.expt standard="Al" eyeball=False calibrated_geom="calibrated.expt"
"""

from sys import exit
from typing import Optional

import numpy as np
import numpy.ma as ma
from matplotlib import pyplot as plt
from matplotlib.widgets import Button, Slider

from iotbx.phil import parse


def module_exists(module_name):
    try:
        __import__(module_name)
    except ModuleNotFoundError as err:
        print(
            err,
            f"""
            This script requires the {module_name} library. Try:
                conda install -c conda-forge {module_name}
            """,
        )
        exit()
    else:
        return True


if module_exists("pyFAI"):
    from pyFAI.calibrant import get_calibrant as pfCalibrant
    from pyFAI.detectors import Detector as pfDetector
    from pyFAI.geometry import Geometry as pfGeometry
    from pyFAI.goniometer import SingleGeometry as pfSingleGeometry
    from pyFAI.gui import jupyter as pfjupyter


# dxtbx and dials must be imported after pyfai,
# the alternative causes pyfai to segv
from dials.util.options import OptionParser, flatten_experiments

phil_scope = parse(
    """
  starting_geom = None
    .type = str
    .help = file containing starting geometry
  standard = None
    .type = str
    .help = calibrant name
  eyeball = True
    .type = bool
    .help = use widget to eyeball a better starting geometry?
  calibrated_geom = None
    .type = str
    .help = file to which the calibrated geometry would be saved
    """
)


def _convert_units(val, unit_in, unit_out):
    """Keep units sanity for pyFAI <--> DIALS conversions
    :param val: the value to be converted
    :param unit_in: Units of input value
    :param unit_out: Units wanted for the input value
    :return: value in new units
    """
    SI = {"A": 1e-10, "nm": 1e-9, "micron": 1e-6, "mm": 1e-3, "cm": 1e-2, "m": 1}
    return val * SI[unit_in] / SI[unit_out]


def parse_command_line():
    """
    Parse command line arguments and read experiments
    """
    usage = (
        "$ dials.powder_calibrate EXPERIMENTS standard_name=standard_name eyeball=True "
        "                                     start_geom_file=file_name final_geom_file [options]"
    )

    parser = OptionParser(
        usage=usage,
        phil=phil_scope,
        check_format=True,
        read_experiments=True,
    )

    parameters, options = parser.parse_args(show_diff_phil=False)
    return parameters, options


def show_fit(geometry, label=None):
    pfjupyter.display(sg=geometry, label=label)
    plt.show()


class Detector(pfDetector):
    def __init__(self, dials_detector):
        px, py = dials_detector.get_pixel_size()
        size_x, size_y = dials_detector.get_image_size()

        super().__init__(
            pixel1=_convert_units(px, "mm", "m"),
            pixel2=_convert_units(py, "mm", "m"),
            max_shape=(size_x, size_y),
        )


class DialsParams:
    def __init__(
        self,
        phil_params,
    ):
        """
        A silly wrapper to make access to useful parameters easier.
        """
        self.start_file = phil_params.input.experiments[0].filename
        self.expt = flatten_experiments(phil_params.input.experiments)[0]
        self.detector = self.expt.detector[0]
        self.beam = self.expt.beam
        self.wavelength = self.beam.get_wavelength()

        self.img_size = self.detector.get_image_size()
        self.image = np.array(self.expt.imageset.get_corrected_data(0)[0]).reshape(
            self.img_size
        )


class Geometry(pfGeometry):
    """
    pyFAI uses the normal from the sample to the detector as the direction
    along which the detector distance and its intersection with the detector
    denotes the detector geometry. This intersection is labelled as `poni` and
    must be differentiated from the intersection of the beam with the detector.
    These two are coincident only in the case of flat (zero tilt) detector.

    Passing geometry back and from DIALS involves passing parameters through
    pyFAI.geometry.getFit2D and PyFAI.geometry.setFit2D which understands
    beam center as a geometry parameter and updates correspondingly the
    pyFAI.geometry and pyFAI.detector objects
    """

    def __init__(self, dials_data, geom_file=None):
        self.dials_data = dials_data
        detector = dials_data.detector

        super().__init__(
            detector=Detector(detector),
            wavelength=_convert_units(dials_data.wavelength, "A", "m"),
        )
        if geom_file:
            # ToDo: but make it dials.expt file

            # read geometry from file
            self.set_param(self.read(geom_file).param)
            beam = self.getFit2D()

            self.beam_x_m = beam["centerX"] / beam["pixelX"]
            self.beam_y_m = beam["centerY"] / beam["pixelY"]

            self.beam_x_px = beam["centerX"]
            self.beam_y_px = beam["centerY"]

        else:
            # fill in geometry from dials_data
            s0 = dials_data.beam.get_s0()

            # beam parameters from dials in mm
            beam_x, beam_y = detector.get_beam_centre(s0)

            # convert to m
            self.beam_x_m = _convert_units(beam_x, "mm", "m")
            self.beam_y_m = _convert_units(beam_y, "mm", "m")

            self.beam_x_px = self.beam_x_m / self.detector.pixel1
            self.beam_y_px = self.beam_y_m / self.detector.pixel2

            self.beam_distance = detector.get_distance()

            # fit2D understands beam parameters
            self.setFit2D(
                directDist=self.beam_distance,
                centerX=self.beam_x_px,
                centerY=self.beam_y_px,
            )

    def update_beam_center(self, beam_coords_px=None, beam_coords_m=None):
        center_x = None
        center_y = None

        if beam_coords_px:
            center_x = beam_coords_px[0]
            center_y = beam_coords_px[1]
        elif beam_coords_m:
            center_x = beam_coords_px[0] / self.detector.pixel1
            center_y = beam_coords_px[1] / self.detector.pixel2

        f2d = self.getFit2D()
        self.setFit2D(directDist=f2d["directDist"], centerX=center_x, centerY=center_y)

    def update_from_ai(self, ai):
        """
        Update geometry with geometry from refined ai
        Is there a more pyFAI friendly of doing this?
        Parameters
        ----------
        ai: ai object
        """
        self.set_param(ai.param)

    def to_phil(self):
        """
        Translate parameters to phil like object such that dials.command_line.modify_geometry can understand them
        """
        phil = [
            "fast_slow_beam_centre="
            + str(int(self.beam_x_px))
            + ","
            + str(int(self.beam_y_px)),
            "distance=" + str(self.beam_distance),
            "wavelength=" + str(self.wavelength),
        ]
        print(phil)
        return phil

    def save_to_expt(self, output="modified.expt"):
        """
        Update the geometry from start_geometry.expt and save to new output
        Pretend DIALS has a python API
        """
        from dials.command_line import modify_geometry

        new_phil = self.to_phil()
        args = [self.dials_data.start_file] + new_phil + ["output=" + output]
        print(args)
        modify_geometry.run(args)

    def __deepcopy__(self, memo=None):
        new = self.__class__(self.dials_data)
        return new


class PowderCalibrator:
    def __init__(
        self,
        dials_data: DialsParams,
        standard: str,
        eyeball: Optional[bool] = True,
        output_geom: Optional[str] = None,
    ):
        """
        Perform geometry calibration using an electron powder standard. Because electron powder
        rings are more dense than X-ray ones, pyFAI struggles to automatically find the correct rings.
        The current method tries to ensure pyFAI has a good starting guess by asking the user to
        eyeball the beam position using a matplotlib widget which overlaps the theoretical calibrator rings
        over the measured ones. The widget part can be turned off setting `eyeball=False`.

        Parameters
        ----------
        dials_data: dials data object
        standard: calibrating standard as str
        eyeball: optional, default=True, toggles the usage of the eyeballing widget.
                When set to False the calibration process only wraps pyFAI and if no
                good starting geometry is given in geom_file then it will probably
                return poor results.
        output_geom: optional, if given it saves the calibrated geometry to file
        """

        self.data = dials_data
        self.eyeball = eyeball
        self.final_geom_file = output_geom

        self.geometry = Geometry(dials_data)

        self.detector = Detector(dials_data.detector)

        self.standard_name = standard
        self.calibrant = pfCalibrant(standard)

        self.calibrant.wavelength = _convert_units(self.data.wavelength, "A", "m")

        self.print_hints()

    def print_hints(self):
        # Tell me which geometry I'm starting from
        print(
            f"Initial geometry from {self.data.start_file}:\n -----\n {self.geometry} \n"
        )

        # Tell me what calibrant I am using
        print(
            f"Starting calibration using {self.standard_name}:\n----- \n {self.calibrant} \n\n "
        )

        # Tell me if I'm using the eyeball widget and how to use
        if self.eyeball:
            print(
                "Using the eyeballing widget.\n"
                "Drag sliders to roughly overlap rings and then save for further fitting. \n"
            )

    def rough_fit_widget(self):
        """
        Matplotlib widget to eyeball the beam center by comparing the theoretical
        and experimental standard.
        """
        fig, ax = plt.subplots()
        ax = pfjupyter.display(
            self.data.image, label="Calibrant overlay on experimental image", ax=ax
        )

        rings_img = self.rings_beam_image(self.geometry, ax)

        ax.set_xlabel("x position [pixels]")
        ax.set_ylabel("y position [pixels]")

        plt.subplots_adjust(left=0.25, bottom=0.25)

        # Make a horizontal slider to control the beam x position.
        ax_beam_x = plt.axes([0.25, 0.1, 0.65, 0.03])

        beam_x_slider = Slider(
            ax=ax_beam_x,
            label="Beam center x [px]",
            valmin=self.detector.max_shape[0] * 0.25,
            valmax=self.detector.max_shape[0] * 0.75,
            valinit=self.geometry.beam_x_px,
        )
        # Make a vertically oriented slider to the beam y position
        ax_beam_y = plt.axes([0.1, 0.25, 0.0225, 0.63])
        beam_y_slider = Slider(
            ax=ax_beam_y,
            label="Beam center y [px]",
            valmin=self.detector.max_shape[1] * 0.25,
            valmax=self.detector.max_shape[1] * 0.75,
            valinit=self.geometry.beam_y_px,
            orientation="vertical",
        )

        def update(val):
            new_geometry = self.geometry.__deepcopy__()

            new_geometry.update_beam_center(
                beam_coords_px=[beam_x_slider.val, beam_y_slider.val]
            )
            rings_img.set_array(self.cal_rings(new_geometry))
            fig.canvas.draw_idle()

        # register the update function with each slider
        beam_x_slider.on_changed(update)
        beam_y_slider.on_changed(update)

        reset_ax = plt.axes([0.8, 0.026, 0.1, 0.04])
        button_res = Button(reset_ax, "Reset", hovercolor="0.975")

        def reset(event):
            beam_x_slider.reset()
            beam_y_slider.reset()

        button_res.on_clicked(reset)

        save_ax = plt.axes([0.5, 0.026, 0.23, 0.04])
        button_save = Button(save_ax, "Save beam and exit", hovercolor="0.975")

        def save_and_exit(event):
            self.geometry.update_beam_center(
                beam_coords_px=[beam_x_slider.val, beam_y_slider.val]
            )
            self.geometry.save_to_expt("eyeballed.expt")

            plt.close()

        button_save.on_clicked(save_and_exit)

        plt.show()

    def rings_beam_image(self, geometry, ax):
        cal_img_masked = self.cal_rings(geometry)

        rings_img = ax.imshow(cal_img_masked, alpha=0.3, cmap="inferno", origin="lower")
        ax.plot(*self.beam_position(geometry), "rx")

        return rings_img

    def beam_position(self, geometry=None):
        """
        Compute beam position
        """
        if geometry:
            beam = [
                geometry.poni2 / self.detector.pixel2,
                geometry.poni1 / self.detector.pixel1,
            ]
        else:
            beam = [
                self.geometry.poni2 / self.detector.pixel2,
                self.geometry.poni1 / self.detector.pixel1,
            ]
        return beam

    def cal_rings(self, geometry=None):
        # for smaller wavelengths (electrons) reduce the rings blurring

        if self.data.wavelength < 1e-1:
            w = 1e-7
        else:
            w = 1e-3

        if geometry:
            cal_img = self.calibrant.fake_calibration_image(geometry, W=w)
        else:
            cal_img = self.calibrant.fake_calibration_image(self.geometry, W=w)

        cal_img_masked = ma.masked_where(cal_img <= 1e-4, cal_img)

        return cal_img_masked

    def calibrate_with_calibrant(
        self,
        num_rings=4,
        fix=("rot1", "rot2", "rot3", "wavelength"),
        verbose=False,
    ):

        if self.eyeball:
            # first use user eyes for rough fit
            self.rough_fit_widget()

        else:
            print(
                "Warning: If the starting geometry is significantly off the fit might return poor results."
            )

        # then use pyFAI for fine calibration
        gonio_geom = pfSingleGeometry(
            label=self.standard_name + " calibrant",
            image=self.data.image,
            calibrant=self.calibrant,
            geometry=self.geometry,
        )

        gonio_geom.extract_cp(max_rings=num_rings)

        if verbose:
            show_fit(gonio_geom, label="Starting geometry")

        gonio_geom.geometry_refinement.refine2(fix=fix)
        ai = gonio_geom.get_ai()

        if verbose:
            show_fit(gonio_geom, label="After pyFAI fit")
            print("Geometry fitted by pyFAI:\n----- \n", ai, "\n")

            # show the cake plot as well
            int2 = ai.integrate2d_ng(
                self.data.image,
                500,
                360,
                unit="q_nm^-1",
            )
            pfjupyter.plot2d(int2, calibrant=self.calibrant)
            plt.title("Are those lines straight?")
            plt.show()

        if self.final_geom_file:
            ai.save(self.final_geom_file)
            print(f"Calibrated geometry saved to {self.final_geom_file}")


if __name__ == "__main__":

    params, opt = parse_command_line()

    standard_params = DialsParams(phil_params=params)

    calibrator = PowderCalibrator(
        dials_data=standard_params,
        standard=params.standard,
        eyeball=params.eyeball,
        output_geom=params.calibrated_geom,
    )

    calibrator.calibrate_with_calibrant(verbose=False)
