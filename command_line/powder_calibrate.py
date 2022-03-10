"""
Calibrate electron powder geometry using X-ray powder diffraction tools.
`pyFAI` is a well established X-ray powder diffraction tool.
See more here: https://doi.org/10.1107/S1600576715004306

The calibration is done in two steps:
    Step 1) A coarse beam position calibration done by the user using a matplolib widget.
    Step 2) Starting from the coarse calibration, pyFAI geometry calibration provides a fine full geometry calibration.
"""

from __future__ import annotations

import logging
import os
from sys import exit
from typing import List, NamedTuple, Optional, Tuple, Union

import numpy as np
import numpy.ma as ma
from matplotlib import pyplot as plt
from matplotlib.widgets import Button, Slider

from libtbx.phil import parse, scope

try:
    import pyFAI  # noqa: F401
except ModuleNotFoundError:
    exit(
        """
            This script requires the pyFAI library. Try:
                $ conda install -c conda-forge pyFAI-base
                   OR if you prefer mamba
                $ mamba install pyFAI-base
            """
    )
from pyFAI.calibrant import get_calibrant as pfCalibrant
from pyFAI.detectors import Detector as pfDetector
from pyFAI.geometry import Geometry as pfGeometry
from pyFAI.goniometer import SingleGeometry as pfSingleGeometry
from pyFAI.gui import jupyter as pfjupyter

# dxtbx and dials must be imported after pyfai,
# the alternative causes pyfai to segv
import dials.util
import dials.util.log
from dials.command_line import modify_geometry
from dials.util.options import ArgumentParser, flatten_experiments
from dials.util.version import dials_version

# Define a logger.
logger = logging.getLogger("dials.command_line.powder_calibrate")

phil_scope = parse(
    """
    output {
        coarse_geom = coarse_geom.expt
            .type = str
            .help = file to which the coarse geometry would be saved
        calibrated_geom = calibrated.expt
            .type = str
            .help = file to which the calibrated geometry would be saved
        pyfai_improvement = pyfai_imporovement.png
            .type = str
            .help = file name for the pyfai calibration effect figure
        straight_lines = straight_lines.png
            .type = str
            .help = file name for the cake plot showing hopefully straight lines
        log = dials.command_name.log
           .type = path
    }
  starting_geom = None
    .type = str
    .help = file containing starting geometry
  standard = None
    .type = str
    .help = calibrant name
  eyeball = True
    .type = bool
    .help = use widget to eyeball a better starting geometry?
    """
)


class Point(NamedTuple):
    """
    Use the blessing of namedtuple for x y coordinates
    """

    slow: float
    fast: float


class ExptParams(NamedTuple):
    """
    Store the .expt parameters that will be used in calibration
    """

    input_file: str
    s0: Tuple[float, float, float]
    beam_on_detector: Point
    wavelength: float
    distance: float
    img_size: Tuple[int, int]
    pix_size: Tuple[float, float]
    image: np.ndarray


class UserArgs(NamedTuple):
    """
    Store the arguments given
    """

    eyeball: bool
    standard: str
    coarse_geom: str
    calibrated_geom: str
    pyfai_improvement: str
    straight_lines: str


def _convert_units(
    val: Union[float, Point], unit_in: str, unit_out: str
) -> Union[float, Point]:
    """Keep units sanity for pyFAI <--> DIALS conversions
    Most parameters will be kept in SI units internally

    :param val: the value/values to be converted
    :param unit_in: Units of input value
    :param unit_out: Units wanted for the input value
    :return: value in new units
    """
    si = {"A": 1e-10, "nm": 1e-9, "micron": 1e-6, "mm": 1e-3, "cm": 1e-2, "m": 1}
    if isinstance(val, Point):
        return Point(*np.array(val) * si[unit_in] / si[unit_out])
    else:
        return val * si[unit_in] / si[unit_out]


def parse_to_tuples(
    args: List[str] = None, phil: scope = phil_scope
) -> Tuple[ExptParams, UserArgs]:
    """
    Parse arguments from command line or not and return named tuple with the useful parameters
    """
    usage = "$ dials.powder_calibrate imported.expt standard=standard_name [options]"

    parser = ArgumentParser(
        usage=usage, phil=phil, read_experiments=True, check_format=True, epilog=__doc__
    )

    params, options = parser.parse_args(args=args, show_diff_phil=False)

    # Configure the logging.
    dials.util.log.config(options.verbose, logfile=params.output.log)

    # Log the dials version
    logger.info(dials_version())

    # Log the difference between the PHIL scope definition and the active PHIL scope,
    # which will include the parsed user inputs.
    diff_phil = parser.diff_phil.as_str()
    if diff_phil:
        logger.info("The following parameters have been modified:\n%s", diff_phil)

    # populate experimental parameters
    experiments = flatten_experiments(params.input.experiments)

    # Complain if command line input is wrong
    if len(experiments) == 0:
        parser.print_help()
        exit("Experiment list required.")
    if len(params.standard) == 0:
        parser.print_help()
        exit("Standard name required.")

    expt = experiments[0]
    detector = expt.detector[0]
    beam = expt.beam

    s0 = beam.get_s0()
    beam_on_detector = detector.get_beam_centre(s0)

    distance = detector.get_distance()
    wavelength = beam.get_wavelength()
    pix_size = detector.get_pixel_size()
    img_size = detector.get_image_size()
    image = np.array(expt.imageset.get_corrected_data(0)[0]).reshape(img_size)

    expt_params = ExptParams(
        input_file=params.input.experiments[0].filename,
        s0=s0,
        beam_on_detector=Point(*beam_on_detector),
        wavelength=wavelength,
        distance=distance,
        img_size=img_size,
        pix_size=pix_size,
        image=image,
    )

    # populate user arguments
    user_args = UserArgs(
        standard=params.standard,
        eyeball=params.eyeball,
        coarse_geom=params.output.coarse_geom,
        calibrated_geom=params.output.calibrated_geom,
        pyfai_improvement=params.output.pyfai_improvement,
        straight_lines=params.output.straight_lines,
    )

    return expt_params, user_args


class Detector(pfDetector):
    def __init__(self, expt_params: ExptParams):
        px, py = expt_params.pix_size
        size_x, size_y = expt_params.img_size

        super().__init__(
            pixel1=_convert_units(px, "mm", "m"),
            pixel2=_convert_units(py, "mm", "m"),
            max_shape=(size_x, size_y),
        )


class Geometry(pfGeometry):
    """
    pyFAI uses the normal from the sample to the detector as the direction
    along which the detector distance and its intersection with the detector
    denotes the detector geometry. This intersection is labelled as `poni` and
    must be differentiated from the intersection of the beam with the detector.
    These two are coincident only in the case of flat (zero tilt) detector.

    Passing geometry back and from DIALS involves passing parameters through
    pyFAI.geometry.setFit2D which understands beam center as a geometry parameter
    and updates correspondingly the pyFAI.geometry and pyFAI.detector objects.
    Similarly, PyFAI.geometry.setFit2D does the reverse.
    """

    def __init__(self, expt_params: ExptParams):
        self.expt_params = expt_params

        super().__init__(
            detector=Detector(expt_params),
            wavelength=_convert_units(expt_params.wavelength, "A", "m"),
        )
        # convert beam parameters from mm to m
        self.beam_m = _convert_units(expt_params.beam_on_detector, "mm", "m")
        self.beam_px = self._beam_to_px()

        self.beam_distance = expt_params.distance

        # pyFAI calibration will need its poni parameters
        self._set_poni()

    def _set_poni(self):
        """
        Call fit2D to translate beam parameters to poni parameters
        """
        self.setFit2D(
            directDist=self.beam_distance,
            centerX=self.beam_px.slow,
            centerY=self.beam_px.fast,
        )

    def _beam_to_px(self):
        """
        Transform beam coordinates from meters to pixels
        """
        if not self.beam_m:
            exit("No beam information provided")
        return Point(
            slow=self.beam_m.slow / self.detector.pixel1,
            fast=self.beam_m.fast / self.detector.pixel2,
        )

    def _beam_to_m(self):
        """
        Transforms beam coordinates from pixels to meters
        """
        if not self.beam_px:
            exit("No beam information provided")
        return Point(
            slow=self.beam_px.slow * self.detector.pixel1,
            fast=self.beam_px.fast * self.detector.pixel2,
        )

    def update_beam_pos(
        self,
        beam_coords_px: Optional[Point] = None,
        beam_coords_m: Optional[Point] = None,
    ):
        if beam_coords_px:
            self.beam_px = beam_coords_px
            self.beam_m = self._beam_to_m()

        elif beam_coords_m:
            self.beam_m = beam_coords_m
            self.beam_px = self._beam_to_px()

        self._set_poni()

    def update_from_ai(self, geom: pfGeometry):
        """
        Update geometry parameters from refined ai

        Parameters
        ----------
        geom: pfGeometry object
        """
        self.set_param(geom.param)
        beam_params = self.getFit2D()

        self.beam_px = Point(slow=beam_params["centerX"], fast=beam_params["centerY"])

        self.beam_m = self._beam_to_m()
        self.beam_distance = beam_params["directDist"]

    def to_parsable(self, only_beam: bool = False) -> List[str]:
        """
        Translate parameters to a parsable list to feed to dials.$command_line_program
        """
        phil = [
            "fast_slow_beam_centre="
            + str(self.beam_px.slow)
            + ","
            + str(self.beam_px.fast)
        ]

        if not only_beam:
            phil += [
                "distance=" + str(self.beam_distance),
                "wavelength=" + str(self.wavelength),
            ]

        return phil

    def save_to_expt(self, output: str | os.PathLike, only_beam: bool = False):
        """
        Update the geometry from start_geometry.expt and save to new output
        Pretend dials.command_line has python API
        """
        new_phil = self.to_parsable(only_beam)
        modify_args = (
            [self.expt_params.input_file] + new_phil + ["output=" + str(output)]
        )
        modify_geometry.run(modify_args)

    def __deepcopy__(self, memo=None):
        new = self.__class__(self.expt_params)
        return new

    def __repr__(self):
        return (
            f"Detector {self.detector.name}: PixelSize= {self.pixel1}, {self.pixel2} m  "
            f"Distance= {self.beam_distance:.2f} mm\n"
            f"Wavelength= {self.wavelength:.3e} m\n"
            f"Beam position on detector: m= {self.beam_m.slow:.4f}m, {self.beam_m.fast:.4f}m  "
            f"px= {self.beam_px.slow:.2f}, {self.beam_px.fast:.2f}"
        )


class EyeballWidget:
    """
    Matplotlib widget to guestimate the beam center by comparing the theoretical
    and experimental standard.
    """

    def __init__(
        self,
        image: np.ndarray,
        start_geometry: Geometry,
        calibrant: pfCalibrant,
        coarse_geom: str,
    ):
        """Eyeball widget
        :param image: 2D array of the image
        :param start_geometry: initial geometry
        :param calibrant: the standard used for calibration
        :param coarse_geom: the file name to which the eyeballed geometry is to saved
        """
        self.image = image
        self.geometry = start_geometry
        self.detector = start_geometry.detector
        self.calibrant = calibrant
        self.coarse_geom = coarse_geom
        self.fig, self.ax = self.set_up_figure()
        self.calibrant_image = self.calibrant_rings_image(self.ax)
        self.beam_slow_slider = self.make_slider("slow")
        self.beam_fast_slider = self.make_slider("fast")

    def __repr__(self):
        calibrant = os.path.splitext(os.path.basename(self.calibrant.filename))[0]
        return f"Eyeballing Widget using {calibrant} Calibrant starting from Geometry: \n {self.geometry}"

    def set_up_figure(self):
        """
        Add title and label axes
        """
        fig, ax = plt.subplots()
        ax = pfjupyter.display(
            self.image,
            label="Calibrant overlay on experimental image",
            ax=ax,
        )
        ax.set_xlabel("slow position [pixels]")
        ax.set_ylabel("fast position [pixels]")
        fig.set_size_inches(10, 10)
        return fig, ax

    def update(self, val):
        """
        Update geometry from slider value
        """
        new_geometry = self.geometry.__deepcopy__()
        new_geometry.update_beam_pos(
            beam_coords_px=Point(self.beam_slow_slider.val, self.beam_fast_slider.val)
        )

        self.calibrant_image.set_array(self.calibrant_rings(new_geometry))
        self.fig.canvas.draw_idle()

    def reset(self, event):
        """
        Reset calibrant image to starting position
        """
        self.beam_slow_slider.reset()
        self.beam_fast_slider.reset()

    def save_and_exit(self, event):
        """
        Save geometry from widget and save to file
        """
        self.geometry.update_beam_pos(
            beam_coords_px=Point(self.beam_slow_slider.val, self.beam_fast_slider.val)
        )
        self.geometry.save_to_expt(
            only_beam=True,
            output=self.coarse_geom,
        )
        plt.close()

    def eyeball(self):
        """
        Update geometry such that beam position is now the center of the moved calibrant rings
        """
        # register the update function with each slider
        self.beam_slow_slider.on_changed(self.update)
        self.beam_fast_slider.on_changed(self.update)

        # register the reset function with reset button
        reset_ax = plt.axes([0.8, 0.026, 0.1, 0.04])
        button_res = Button(reset_ax, "Reset", hovercolor="0.975")
        button_res.on_clicked(self.reset)

        # register the save and exit function with save button
        save_ax = plt.axes([0.5, 0.026, 0.23, 0.04])
        button_save = Button(save_ax, "Save beam and exit", hovercolor="0.975")
        button_save.on_clicked(self.save_and_exit)

        # finally, show widget
        plt.show()

    def make_slider(self, label):
        slider = None
        if label == "slow":
            plt.subplots_adjust(bottom=0.25)
            # Make a horizontal slider to control the beam slow position.
            ax_beam_slow = plt.axes([0.25, 0.1, 0.65, 0.03])

            slider = Slider(
                ax=ax_beam_slow,
                label="Beam center slow [px]",
                valmin=self.detector.max_shape[0] * 0.25,
                valmax=self.detector.max_shape[0] * 0.75,
                valinit=self.geometry.beam_px.slow,
            )

        elif label == "fast":
            plt.subplots_adjust(left=0.25)
            # Make a vertically oriented slider to the beam fast position
            ax_beam_fast = plt.axes([0.1, 0.25, 0.0225, 0.63])
            slider = Slider(
                ax=ax_beam_fast,
                label="Beam center fast [px]",
                valmin=self.detector.max_shape[1] * 0.25,
                valmax=self.detector.max_shape[1] * 0.75,
                valinit=self.geometry.beam_px.fast,
                orientation="vertical",
            )
        return slider

    def calibrant_rings_image(self, ax: plt.axes) -> plt.axes:
        cal_img_masked = self.calibrant_rings()

        rings_img = ax.imshow(cal_img_masked, alpha=0.3, cmap="inferno", origin="lower")
        # add the beam position
        ax.plot(*self.beam_position(self.geometry), "rx")

        return rings_img

    def beam_position(self, geometry: Optional[Geometry] = None) -> Point:
        """
        Compute the beam position either from a new geometry or from self.geometry
        """
        if geometry:
            beam = Point(
                slow=geometry.poni2 / self.detector.pixel2,
                fast=geometry.poni1 / self.detector.pixel1,
            )
        else:
            beam = Point(
                slow=self.geometry.poni2 / self.detector.pixel2,
                fast=self.geometry.poni1 / self.detector.pixel1,
            )
        return beam

    def calibrant_rings(self, geometry: Optional[Geometry] = None) -> np.ndarray:
        """
        Make a masked array for the calibration rings. If a geometry parameters
        is given then use that geometry for the calibrant, otherwise use the
        class geometry.
        For smaller wavelengths (electrons) reduce the rings blurring
        """

        if self.calibrant.wavelength < 1e-1:
            w = 1e-7
        else:
            w = 1e-3

        if geometry:
            cal_img = self.calibrant.fake_calibration_image(geometry, W=w)
        else:
            cal_img = self.calibrant.fake_calibration_image(self.geometry, W=w)

        cal_img_masked = ma.masked_where(cal_img <= 1e-4, cal_img)

        return cal_img_masked


class PowderCalibrator:
    def __init__(
        self,
        starting_geom: Optional[str] = None,
        standard: Optional[str] = None,
        eyeball: bool = True,
        coarse_geom: str = "coarse_geom.expt",
        calibrated_geom: str = "calibrated.expt",
        pyfai_improvement: str = "pyfai_improvement.png",
        straight_lines: str = "straight_lines.png",
    ):
        """
        Perform geometry calibration using an electron powder standard. Because electron powder
        rings are more dense than X-ray ones, pyFAI struggles to automatically find the correct rings.
        The current method tries to ensure pyFAI has a good starting guess by asking the user to
        eyeball the beam position using a matplotlib widget which overlaps the theoretical calibrator rings
        over the measured ones. The widget part can be turned off setting `eyeball=False`.

        :param starting_geom: str
                expt file containing starting geometry
        :param standard: str
                calibrating standard
        :param eyeball: bool
                optional, default=True, toggles the usage of the eyeballing widget.
                When set to False the calibration process only wraps pyFAI and if no
                good starting geometry is given in geom_file then it will probably
                return poor results.
        :param coarse_geom: str
                optional, if given saves the coarse geometry to this file
        :param calibrated_geom: str
                optional, if given saves the calibrated geometry to this file
        param pyfai_improvement: str
                optional, file name used for saving the before and after pyfai calibration
        param straight_lines: str
                optional, file name used for saving the cake plot, if fir was succesful
                these lines should be straight

        Examples
        --------
        1. Use the widget to generate a starting geometry. Overwrite the default name of the
         `coarse_geom.expt`. Do this from command line:
            $ dials.powder_calibrate poor_geom.expt standard="Al" coarse_geom="coarse_geom_new.expt"

        2. Starting from a coarse geometry file, perform fine calibration with pyFAI bypassing the
        widget tool. Do this using python API:
            >>> from dials.command_line.powder_calibrate import PowderCalibrator
            >>> al_calibrator = PowderCalibrator("coarse_geom.expt", standard="Al", eyeball=False)
            >>> al_calibrator.calibrate_with_calibrant()

        3. Do the same from command line:
            $ dials.powder_calibrate coarse_geom.expt standard="Al" eyeball=False
        """
        if starting_geom and standard:
            # API usage
            args = [
                str(starting_geom),
                "standard=" + standard,
                "eyeball=" + str(eyeball),
                "coarse_geom=" + coarse_geom,
                "calibrated_geom=" + calibrated_geom,
                "pyfai_improvement=" + pyfai_improvement,
                "straight_lines=" + straight_lines,
            ]
            self.expt_params, self.user_args = parse_to_tuples(args=args)
        elif not starting_geom and not standard:
            # command line usage
            self.expt_params, self.user_args = parse_to_tuples()
        elif not standard:
            exit("Standard missing")
        elif not starting_geom:
            exit("Starting geometry missing")

        self.geometry = Geometry(self.expt_params)
        self.detector = self.geometry.detector

        # ToDo: calibrant class?
        self.calibrant = pfCalibrant(self.user_args.standard)
        self.calibrant.wavelength = _convert_units(
            self.expt_params.wavelength, "A", "m"
        )

        self.print_hints()

    def __repr__(self):
        return (
            f"{self.user_args.standard} Calibrator for {self.expt_params.input_file} \n"
            f"Current geometry: \n {self.geometry}"
        )

    def print_hints(self):
        # Tell me which geometry I'm starting from
        logger.info(
            f"Initial geometry from {self.expt_params.input_file}:\n-----\n {self.geometry} \n"
        )

        # Tell me what calibrant I am using
        logger.info(
            f"Starting calibration using {self.user_args.standard}:\n-----\n {self.calibrant} \n "
        )

        # Tell me if I'm using the EyeballWidget and how to use
        if self.user_args.eyeball:
            logger.info(
                "Using the eyeballing widget.\n"
                "Drag sliders to roughly overlap rings and then save for further fitting. \n"
            )

    def show_pyfai_improvement(
        self,
        before_geometry: pfSingleGeometry,
        after_geometry: pfSingleGeometry,
        show: bool = True,
    ):
        fig, (ax1, ax2) = plt.subplots(1, 2)
        fig.set_size_inches(20, 8)
        pfjupyter.display(sg=before_geometry, label="Initial geometry", ax=ax1)
        pfjupyter.display(sg=after_geometry, label="After pyFAI fit", ax=ax2)

        if show:
            plt.show()
        else:
            fig.savefig(self.user_args.pyfai_improvement, bbox_inches="tight")
            plt.close(fig)

    def show_straight_lines(self, ai: pfGeometry, show: bool = True):
        # show the cake plot as well
        int2 = ai.integrate2d_ng(
            data=self.expt_params.image,
            npt_rad=500,
            npt_azim=360,
            unit="q_nm^-1",
        )
        fig, ax = plt.subplots()
        fig.set_size_inches(8, 8)
        pfjupyter.plot2d(int2, calibrant=self.calibrant, ax=ax)
        plt.title("Are those lines straight?")

        if show:
            plt.show()
        else:
            fig.savefig(self.user_args.straight_lines, bbox_inches="tight")
            plt.close(fig)

    def calibrate_with_calibrant(
        self,
        num_rings: int = 4,
        fix: Tuple = ("rot1", "rot2", "rot3", "wavelength"),
        plots: bool = True,
    ):

        if self.user_args.eyeball:
            # first use user eyes for rough fit
            eyeballing_widget = EyeballWidget(
                image=self.expt_params.image,
                start_geometry=self.geometry,
                calibrant=self.calibrant,
                coarse_geom=self.user_args.coarse_geom,
            )
            eyeballing_widget.eyeball()

        else:
            logger.warning(
                "If the starting geometry is significantly off the fit might return poor results."
            )

        # then use pyFAI for fine calibration
        gonio_geom = pfSingleGeometry(
            label=self.user_args.standard + " calibrant",
            image=self.expt_params.image,
            calibrant=self.calibrant,
            geometry=self.geometry,
        )

        gonio_geom.extract_cp(max_rings=num_rings)
        starting_geom = gonio_geom

        gonio_geom.geometry_refinement.refine2(fix=fix)

        # generate plot showing the pyFAI refinement effect
        self.show_pyfai_improvement(
            before_geometry=starting_geom, after_geometry=gonio_geom, show=plots
        )

        # update geometry and save to .expt
        ai = gonio_geom.get_ai()
        self.geometry.update_from_ai(ai)
        self.geometry.save_to_expt(output=self.user_args.calibrated_geom)

        logger.info(f"Geometry fitted by pyFAI:\n----- \n {self.geometry} \n")
        self.show_straight_lines(ai, show=plots)


@dials.util.show_mail_handle_errors()
def run(args=None):
    calibrator = PowderCalibrator()
    calibrator.calibrate_with_calibrant(plots=False)


if __name__ == "__main__":
    run()
