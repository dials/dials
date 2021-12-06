"""
Calibrate geometry using powder tools.
    `pyFAI` is a well established X-ray powder diffraction tool.
https://doi.org/10.1107/S1600576715004306
"""

from sys import exit

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.widgets import Button, Slider

from iotbx.phil import parse


def module_exists(module_name):
    try:
        __import__(module_name)
    except ModuleNotFoundError as err:
        if module_name == "circle_fit":
            module_name = "circle-fit"
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
from dxtbx.model.experiment_list import ExperimentList as experiment_list

from dials.array_family import flex
from dials.util.options import OptionParser, flatten_experiments, flatten_reflections

phil_scope = parse(
    """
  calibrant = None
    .type = str
    .help = calibrant name
  poni = None
     .type = str
     .help = poni file containing geometry
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


class Detector(pfDetector):
    def __init__(self, dials_detector):
        px, py = dials_detector.get_pixel_size()
        size_x, size_y = dials_detector.get_image_size()

        super().__init__(
            pixel1=_convert_units(px, "mm", "m"),
            pixel2=_convert_units(py, "mm", "m"),
            max_shape=(size_x, size_y),
        )


class Geometry(pfGeometry):
    """
    pyFAI uses the normal from the sample to the detector
    as a geometry parameter, labelled as poni instead of the beam center.
    These two are coincident only in the case of flat (zero tilt) detector.

    Passing geometry back and from involves passing parameters through
    pyFAI.geometry.getFit2D and PyFAI.geometry.setFit2D which understands
    beam center as a geometry parameter and updates correspondingly the
    pyFAI.geometry and pyFAI.detector objects
    """

    def __init__(self, dials_data, poni_file=None):
        self.dials_data = dials_data

        detector = dials_data.detector
        super().__init__(
            detector=Detector(detector),
            wavelength=_convert_units(dials_data.wavelength, "A", "m"),
        )

        s0 = dials_data.beam.get_s0()

        # beam parameters from dials in mm
        beam_x, beam_y = detector.get_beam_centre(s0)

        # convert to m
        self.beam_x_m = _convert_units(beam_x, "mm", "m")
        self.beam_y_m = _convert_units(beam_y, "mm", "m")

        self.beam_x_px = self.beam_x_m / self.detector.pixel1
        self.beam_y_px = self.beam_y_m / self.detector.pixel2

        beam_distance = detector.get_distance()

        # fit2D understands beam parameters
        self.setFit2D(
            directDist=beam_distance, centerX=self.beam_x_px, centerY=self.beam_y_px
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

    def __deepcopy__(self, memo=None):
        new = self.__class__(self.dials_data)
        return new


def parse_command_line():
    """
    Parse command line arguments and read experiments
    """
    usage = "$ dials.powder_calibrate EXPERIMENTS REFLECTIONS [options]"

    parser = OptionParser(
        usage=usage,
        phil=phil_scope,
        check_format=True,
        read_reflections=True,
        read_experiments=True,
    )

    params, options = parser.parse_args(show_diff_phil=False)
    return params, options


def show_fit(geometry, label=None):
    pfjupyter.display(sg=geometry, label=label)
    plt.show()


class DialsData:
    def __init__(self, params=None, expt_file=None, refl_file=None, num_images=1):
        self.params = params

        self.expt_file = expt_file
        self.refl_file = refl_file
        self.num_image = num_images

        self.expts = self._expt()
        self.refls = self._refl()

        self.detector = self.expts[0].detector[0]
        self.beam = self.expts[0].beam
        self.wavelength = self.beam.get_wavelength()

        self.img_size = self.detector.get_image_size()
        self.image = np.array(self.expts[0].imageset.get_corrected_data(0)[0]).reshape(
            self.img_size
        )

    def _expt(self):
        experiments = None
        if self.params:
            experiments = flatten_experiments(self.params.input.experiments)
        elif self.expt_file:
            experiments = experiment_list.from_file(self.expt_file)
        else:
            exit("No experiments file was given")
        return experiments

    def _refl(self):
        reflections_table = None
        if self.params:
            reflections_table = flatten_reflections(self.params.input.reflections)
        elif self.refl_file:
            reflections_table = flex.reflection_table.from_file(self.refl_file)
        return reflections_table


class PowderCalibrator:
    def __init__(self, dials_data, standard, poni_file=None):
        self.data = dials_data
        self.poni_file = poni_file

        self.geometry = Geometry(dials_data, poni_file)
        self.detector = Detector(dials_data.detector)

        self.standard_name = standard
        self.calibrant = pfCalibrant(standard)
        self.calibrant.wavelength = _convert_units(self.data.wavelength, "A", "m")

        print("Initial geometry:\n -----\n", self.geometry, "\n")
        print("Start calibration using:\n----- \n", self.calibrant, "\n")

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
            print("Rough estimate geometry:\n ----\n", self.geometry, "\n")
            plt.close()

        button_save.on_clicked(save_and_exit)

        plt.show()

    def rings_beam_image(self, geometry, ax):
        cal_img_masked = self.cal_rings(geometry)

        rings_img = ax.imshow(cal_img_masked, alpha=0.3, cmap="inferno", origin="lower")
        ax.plot(*self.beam_position(geometry), "rx")

        return rings_img

    def beam_position(self, geometry=None):
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
        # for smaller wavelengths reduce the rings blurring
        if self.data.wavelength < 1e-1:
            w = 1e-7
        else:
            w = 1e-3

        if geometry:
            cal_img = self.calibrant.fake_calibration_image(geometry, W=w)
        else:
            cal_img = self.calibrant.fake_calibration_image(self.geometry, W=w)
        cal_img_masked = np.ma.masked_where(cal_img <= 1e-4, cal_img)
        return cal_img_masked

    def calibrate_with_calibrant(
        self,
        num_rings=4,
        fix=("rot1", "rot2", "rot3", "wavelength"),
        ref=False,
        verbose=False,
    ):

        # first use user eyes for rough fit
        self.rough_fit_widget()

        # then use pyFAI for fine calibration
        sg = pfSingleGeometry(
            label=self.standard_name + " calibrant",
            image=self.data.image,
            calibrant=self.calibrant,
            geometry=self.geometry,
        )

        sg.extract_cp(max_rings=num_rings)
        if verbose:
            show_fit(sg, label="Before pyFAI fit")

        if ref:
            sg.geometry_refinement.refine2(fix=fix)
            if verbose:
                show_fit(sg, label="After pyFAI fit")
                print("Geometry fitted by pyFAI:\n----- \n", sg.get_ai(), "\n")

            ai = sg.get_ai()
            int2 = ai.integrate2d_ng(
                self.data.image, 500, 360, unit="q_nm^-1", filename="integrated.edf"
            )
            if verbose:
                pfjupyter.plot2d(int2, calibrant=self.calibrant)
                plt.show()


if __name__ == "__main__":

    Al_data = DialsData(
        expt_file="imported.expt",
    )

    calibrator = PowderCalibrator(Al_data, "Al")
    calibrator.calibrate_with_calibrant(ref=True, verbose=True)
