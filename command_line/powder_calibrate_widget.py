"""
Calibrate geometry using powder tools.
    `pyFAI` is a well established X-ray powder diffraction tool.
https://doi.org/10.1107/S1600576715004306
    `circle-fit` is a nifty module that fits an ellipse to a set of points,
even if the points are covering a section as small as a quarter of the conic.
https://doi.org/10.1016/j.csda.2010.12.012
NOTE: The algorithms assumes that find_spots worked well at separating
signal from noise especially close to the beam, assumption which might be hard
to achieve for electron data.
"""

from sys import exit

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.widgets import Button, Slider

import dxtbx.flumpy as flumpy
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
    from pyFAI.gui import jupyter as pfjupyter


# dxtbx and dials must be imported after pyfai,
# the alternative causes pyfai to segv
from dxtbx.model.experiment_list import ExperimentList as experiment_list

from dials.array_family import flex

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
    """Keep units sanity for pyFAI <--> DIALS conversions.
    :param val: the value to be converted.
    :param unit_in: Units of input value.
    :param unit_out: Units wanted for the input value.
    :return: value in new units.
    """
    SI = {"A": 1e-10, "nm": 1e-9, "micron": 1e-6, "mm": 1e-3, "cm": 1e-2, "m": 1}
    return val * SI[unit_in] / SI[unit_out]


class Dials_pfDetector(pfDetector):
    def __init__(self, dials_detector, name="detector"):
        px, py = dials_detector.get_pixel_size()
        size_x, size_y = dials_detector.get_image_size()

        # horrible hack, should use the factory correctly instead
        Dials_pfDetector.__name__ = name
        super().__init__(
            pixel1=_convert_units(px, "mm", "m"),
            pixel2=_convert_units(py, "mm", "m"),
            max_shape=(size_x, size_y),
        )


class Dials_pfGeometry(pfGeometry):
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
        detector = dials_data.detector
        super().__init__(
            detector=Dials_pfDetector(detector),
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


def parse_command_line():
    """
    Parse command line arguments and read experiments
    """

    from dials.util.options import OptionParser

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


def show_fit(geometry):
    pfjupyter.display(sg=geometry)
    plt.show()


def radial_distance(spots_coords: np.array, beam_coords: np.array) -> np.array:
    """
    The radial distance to spot from beam position.
    Parameters
    ----------
    spots_coords: array, shape=(n, 2)
        position (x,y) of spots in px
    beam_coords : array, shape=(2,)
        beam (x, y) position on detector in px
    Returns
    -------
        rad_distances: array, shape=(n,2)
    """
    # rad_distances = np.sqrt(np.sum((spots_coords - beam_coords) ** 2, axis=1))
    x0, y0 = beam_coords
    x, y = spots_coords.T
    rad_distances = np.hypot(x - x0, y - y0)
    return rad_distances


def peaks_to_spectrum(dials_data, beam):
    coords = flumpy.to_numpy(dials_data.refls["xyzobs.px.value"])

    # bin by integer pixels
    coords_xy = coords[:, :2]
    rad_dis = radial_distance(coords_xy, beam).astype(int)

    intensities = flumpy.to_numpy(dials_data.refls["intensity.sum.value"])
    weighted_hist = np.bincount(rad_dis, weights=intensities)
    unweighted_histogram = np.bincount(rad_dis)
    spectrum = np.divide(
        weighted_hist,
        unweighted_histogram,
        out=np.zeros_like(weighted_hist),
        where=unweighted_histogram != 0,
    )
    return spectrum


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

        self.shoeboxes = self._shoeboxes_coords()

        self.img_size = self.detector.get_image_size()
        self.image = np.array(self.expts[0].imageset.get_corrected_data(0)[0]).reshape(
            self.img_size
        )
        self.sum_max_images = self.stack_images_sum_max()

    def _expt(self):
        if self.params:
            from dials.util.options import flatten_experiments

            experiments_list = flatten_experiments(self.params.input.experiments)
        elif self.expt_file:
            experiments_list = experiment_list.from_file(self.expt_file)
        else:
            exit("No experiments file was given")
        return experiments_list

    def _refl(self):
        if self.params:
            from dials.util.options import flatten_reflections

            reflections_table = flatten_reflections(self.params.input.reflections)
        elif self.refl_file:
            reflections_table = flex.reflection_table.from_file(self.refl_file)
        else:
            exit("No reflection file was given")
        return reflections_table

    def stack_images_sum_max(self):
        """
        Aggregate all images in the dataset by only keeping the max for each pixel
        Returns
        -------
        stacked_images
        """
        stacked_images = np.zeros_like(self.image)
        for i in range(self.num_image):
            image = np.array(self.expts[0].imageset.get_corrected_data(0)[i]).reshape(
                self.image.shape
            )
            # stacked_images = pyFAI.average.average_images(data_list, filter_='max', fformat=None)
            stacked_images = np.maximum(stacked_images, image)
        return stacked_images

    def _shoeboxes_coords(self):

        sbpixels = []
        for i in range(len(self.refls)):
            shoeboxes = self.refls[i]["shoebox"]
            sbpixels.append(flumpy.to_numpy(shoeboxes.coords()))
        return np.vstack(sbpixels)[:, :2]


class PowderCalibrator:
    def __init__(self, dials_data, standard, poni_file=None):
        self.data = dials_data
        self.poni_file = poni_file

        self.geometry = Dials_pfGeometry(dials_data, poni_file)
        print("Initial geometry:\n -----\n", self.geometry, "\n")

        self.detector = Dials_pfDetector(dials_data.detector)

        self.standard_name = standard
        self.calibrant = pfCalibrant(standard)
        self.calibrant.wavelength = _convert_units(self.data.wavelength, "A", "m")
        print("\n Start calibration using:\n----- \n", self.calibrant, "\n")

    def rough_fit_widget(self):
        """
        Matplotlib widget to eyeball the beam center by comparing the theoretical
        and experimental standard.
        """
        fig, ax = plt.subplots()
        ax = pfjupyter.display(
            self.data.image, label="Calibrant overlay on experimental image", ax=ax
        )

        rings = self.rings_image(self.geometry, ax)

        ax.set_xlabel("x position [pixels]")
        ax.set_ylabel("y position [pixels]")

        plt.subplots_adjust(left=0.25, bottom=0.25)

        # Make a horizontal slider to control the beam x position.
        ax_beamx = plt.axes([0.25, 0.1, 0.65, 0.03])
        beamx_slider = Slider(
            ax=ax_beamx,
            label="Beam center x [px]",
            valmin=self.detector.max_shape[0] * 0.25,
            valmax=self.detector.max_shape[0] * 0.75,
            valinit=self.geometry.beam_x_px,
        )
        # Make a vertically oriented slider to the beam y position
        ax_beamy = plt.axes([0.1, 0.25, 0.0225, 0.63])
        beamy_slider = Slider(
            ax=ax_beamy,
            label="Beam center y [px]",
            valmin=self.detector.max_shape[1] * 0.25,
            valmax=self.detector.max_shape[1] * 0.75,
            valinit=self.geometry.beam_y_px,
            orientation="vertical",
        )

        def update(val):
            self.geometry.update_beam_center(
                beam_coords_px=[beamx_slider.val, beamy_slider.val]
            )
            rings.set_array(self.cal_rings(self.geometry))
            fig.canvas.draw_idle()

        # register the update function with each slider
        beamx_slider.on_changed(update)
        beamy_slider.on_changed(update)

        reset_ax = plt.axes([0.8, 0.025, 0.1, 0.04])
        button_res = Button(reset_ax, "Reset", hovercolor="0.975")

        def reset(event):
            beamx_slider.reset()
            beamy_slider.reset()

        button_res.on_clicked(reset)

        save_ax = plt.axes([0.5, 0.025, 0.23, 0.04])
        button_save = Button(save_ax, "Save beam and exit", hovercolor="0.975")

        def save_and_exit(event):
            print("Saved geometry:\n ----\n", self.geometry, "\n")
            plt.close()

        button_save.on_clicked(save_and_exit)

        plt.show()

    def update_geometry(self, beam_px=None):
        # update geometry using the Fit2D route
        f2d = self.geometry.getFit2D()
        self.geometry.setFit2D(
            directDist=f2d["directDist"], centerX=beam_px[0], centerY=beam_px[1]
        )
        print("Updated model geometry to", self.geometry)

    def rings_image(self, geometry, ax):
        if self.data.wavelength < 1e-1:
            w = 1e-7
        else:
            w = 1e-3

        cal_img = self.calibrant.fake_calibration_image(geometry, W=w)
        cal_img_masked = np.ma.masked_where(cal_img <= 1e-4, cal_img)

        rings = ax.imshow(cal_img_masked, alpha=0.3, cmap="inferno", origin="lower")

        beam = [
            self.geometry.poni2 / self.detector.pixel2,
            self.geometry.poni1 / self.detector.pixel1,
        ]
        ax.plot(*beam, "rx")

        return rings

    def cal_rings(self, geometry):
        # for smaller wavelengths reduce the blurring of the rings
        if self.data.wavelength < 1e-1:
            w = 1e-7
        else:
            w = 1e-3

        cal_img = self.calibrant.fake_calibration_image(geometry, W=w)
        cal_img_masked = np.ma.masked_where(cal_img <= 1e-4, cal_img)
        return cal_img_masked


if __name__ == "__main__":

    Al_data = DialsData(
        expt_file="/home/fyi77748/Data/Al_standard/imported.expt",
        refl_file="/home/fyi77748/Data/Al_standard/strong.refl",
    )

    # params, opt = parse_command_line()
    # Al_data = DialsData(params=params)

    calibrator = PowderCalibrator(Al_data, "Al")
    calibrator.rough_fit_widget()
