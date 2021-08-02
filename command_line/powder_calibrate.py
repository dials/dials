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
from matplotlib.pyplot import subplots
from scipy.signal import find_peaks

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
    from pyFAI.azimuthalIntegrator import AzimuthalIntegrator as pfAI
    from pyFAI.calibrant import get_calibrant as pfCalibrant
    from pyFAI.detectors import Detector as pfDetector
    from pyFAI.geometry import Geometry as pfGeometry
    from pyFAI.goniometer import SingleGeometry as pfSingleGeometry
    from pyFAI.gui import jupyter as pfjupyter

if module_exists("circle_fit"):
    import circle_fit as cf

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
    :return: Value in new units.
    """
    SI = {"A": 1e-10, "nm": 1e-9, "microm": 1e-6, "mm": 1e-3, "cm": 1e-2, "m": 1}
    return val * SI[unit_in] / SI[unit_out]


class Dials_pfDetector(pfDetector):
    def __init__(self, dials_detector):
        px, py = dials_detector.get_pixel_size()
        size_x, size_y = dials_detector.get_image_size()

        # horrible hack, should use the factory correctly instead
        Dials_pfDetector.__name__ = "detector"
        super().__init__(
            pixel1=_convert_units(px, "mm", "m"),
            pixel2=_convert_units(py, "mm", "m"),
            max_shape=(size_x, size_y),
        )


# def _make_pfDetector(detector):
#     """
#     make a PyFAI Detector object from DIALS experiment
#     """
#     px, py = detector.get_pixel_size()
#     size_x, size_y = detector.get_image_size()
#
#     return pfDetector(
#         pixel1=_convert_units(px, "mm", "m"),
#         pixel2=_convert_units(py, "mm", "m"),
#         max_shape=(size_x, size_y),
#     )


class Dials_pfGeometry(pfGeometry):
    def __init__(self, dials_data, poni_file=None):
        detector = dials_data.detector
        s0 = dials_data.beam.get_s0()
        x, y = detector.get_beam_centre(s0)
        distance = detector.get_distance()

        super().__init__(
            dist=_convert_units(distance, "mm", "m"),
            poni1=_convert_units(y, "mm", "m"),
            poni2=_convert_units(x, "mm", "m"),
            detector=Dials_pfDetector(detector),
            wavelength=_convert_units(dials_data.wavelength, "A", "m"),
        )


# def _make_pfGeometry(detector, s0, wavelength, poni_file=None):
#     """
#     make a PyFAI geometry object from DIALS experiment
#     """
#     if poni_file:
#         pass
#     else:
#         x, y = detector.get_beam_centre(s0)
#         distance = detector.get_distance()
#         geom = pfGeometry(
#             dist=_convert_units(distance, "mm", "m"),
#             poni1=_convert_units(y, "mm", "m"),
#             poni2=_convert_units(x, "mm", "m"),
#             detector=_make_pfDetector(detector),
#             wavelength=wavelength,
#         )
#     return geom


def _make_pfAI(detector, geometry, wavelength):
    return pfAI(
        dist=geometry.dist,
        poni1=geometry.poni1,
        poni2=geometry.poni2,
        detector=detector,
        wavelength=wavelength,
    )


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


def radial_distance(spots_coords, beam_coords):
    """
    The radial distance to spot from beam position.

    Parameters
    ----------
    spots_coords: ndarray, shape=(n, 2)
        position (x,y) of spots in px
    beam_coords : ndarray, shape=(2,)
        beam (x, y) position on detector in px

    Returns
    -------
        rad_distances: ndarray, shape=(n,2)
    """

    rad_distances = np.sqrt(np.sum((spots_coords - beam_coords) ** 2, axis=1))
    return rad_distances


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
    def __init__(self, dials_data, poni_file=None):
        self.data = dials_data
        self.poni_file = poni_file

        self.geometry = Dials_pfGeometry(dials_data, poni_file)
        print("Initial geometry:\n -----\n", self.geometry, "\n")

        self.detector = Dials_pfDetector(dials_data.detector)

    def beam_from_fit_circle(self, rings="all", shoebox=False):
        """
        Find beam position by fitting a conic to data
        Parameters
        ----------
        rings: str
            options:'all', 'first'
            all: fit one conic to all data for an initial guess
            first: use only spots at high resolution,
                   ideally these are scattered around a conic
            default: 'all'
        shoebox: bool
            'True': use shoeboxes pixels as data to fit to
            'False': use instead centroids
            default: False
        Returns
        -------
        beam_coords
        """

        if shoebox:
            coords = self.shoebox_coords()
            coords_xy = coords[:, :2]
        else:
            coords = flumpy.to_numpy(self.data.refls["xyzobs.px.value"])
            coords_xy = coords[:, :2]

        if rings == "all":
            xc_f, yc_f, r, s = cf.hyper_fit(coords_xy)
            beam = [xc_f, yc_f]
            print("Beam position on detector from circle fit:", beam)

            fig = plt.figure()
            circle = plt.Circle((xc_f, yc_f), r, color="b", fill=False)
            ax = fig.gca()
            ax.add_patch(circle)

            plt.plot(*beam, "bx")
            plt.imshow(self.data.image)
            plt.scatter(*zip(*coords_xy))

        elif rings == "first":
            peaks = self.first_circle_lims()
            pix_rad = peaks[1] / self.detector.pixel1
            beam = self.beam
            first_circle = radial_distance(coords_xy, beam) <= (pix_rad - 5)
            xc_f, yc_f, r, s = cf.hyper_fit(coords_xy[first_circle])
            beam = [xc_f, yc_f]
            self.beam = beam
            print("Beam position on detector from circle fit:", beam)

            fig = plt.figure()
            circle = plt.Circle((xc_f, yc_f), r, color="b", fill=False)
            ax = fig.gca()
            ax.add_patch(circle)

            plt.plot(*beam, "bx")
            plt.imshow(self.data.image)
            plt.scatter(*zip(*coords_xy[first_circle]))

        else:
            exit(f"{rings} value not understood.")
        plt.show()

        self.beam = beam
        return beam

    def first_circle_lims(self):
        ai = _make_pfAI(self.detector, self.geometry, self.data.wavelength)

        int = ai.integrate1d(self.data.image, self.detector.max_shape[0], unit="r_mm")
        pfjupyter.plot1d(int)
        peaks = self.find_peaks(int)
        print("Found peaks:", peaks)
        plt.scatter(int.radial[peaks], int.intensity[peaks])
        plt.show()
        return _convert_units(int.radial[peaks], "mm", "m")

    def find_peaks(self, int):
        peaks, _ = find_peaks(int.intensity, prominence=3, width=[1, 50])
        return peaks

    def update_geometry(self, beam_px=None):
        if beam_px:
            self.geometry.poni1 = beam_px[1] * self.detector.pixel1
            self.geometry.poni2 = beam_px[0] * self.detector.pixel2
        print(self.geometry)

    def set_calibrant(self, name):
        self.calibrant_name = name
        self.calibrant = pfCalibrant(name)
        self.calibrant.wavelength = _convert_units(self.data.wavelength, "A", "m")
        print("\n Using:\n----- \n", self.calibrant)

    def calibrate_with_calibrant(
        self,
        name,
        num_rings=4,
        fix=["rot1", "rot2", "rot3", "wavelength"],
        ref=False,
        verbose=False,
    ):
        self.set_calibrant(name)
        if verbose:
            self.show_calibrant(self.geometry)

        sg = pfSingleGeometry(
            label=name + " calibrant",
            image=self.data.image,
            calibrant=self.calibrant,
            geometry=self.geometry,
        )

        sg.extract_cp(max_rings=num_rings)
        if verbose:
            show_fit(sg)

        if ref:
            sg.geometry_refinement.refine2(fix=fix)
            if verbose:
                show_fit(sg)

            ai = sg.get_ai()
            int2 = ai.integrate2d_ng(
                self.data.image, 500, 360, unit="q_nm^-1", filename="integrated.edf"
            )
            if verbose:
                pfjupyter.plot2d(int2, calibrant=self.calibrant)
                plt.show()

    def tricky_fit_with_calibrant(
        self,
        name,
        num_rings=4,
        fix=["rot1", "rot2", "rot3", "wavelength"],
        verbose=False,
    ):
        self.calibrate_with_calibrant(name=name, num_rings=4, fix=fix, verbose=verbose)
        beam = self.beam_from_fit_circle()
        self.update_geometry(beam_px=beam)
        self.calibrate_with_calibrant(
            name=name, num_rings=4, fix=fix, ref=True, verbose=verbose
        )

    def show_calibrant(self, geometry):
        fig, ax = subplots(1, 2, figsize=(8, 5))
        pfjupyter.display(self.data.image, label="Experiment image", ax=ax[0])

        # for smaller wavelengths reduce the blurring of the rings
        if self.data.wavelength < 1e-1:
            w = 1e-7
        else:
            w = 1e-3

        cal_img = self.calibrant.fake_calibration_image(geometry, W=w)
        pfjupyter.display(cal_img, label=self.calibrant_name + " calibrant", ax=ax[1])
        plt.show()

    def show_fit(self, geometry):
        pfjupyter.display(sg=geometry)
        plt.show()

    def save_refined_geometry(self, file_name):
        if self.sg:
            self.sg.geometry_refinement.save(file_name)
        else:
            pass


if __name__ == "__main__":

    Al_data = DialsData(
        expt_file="/home/elena/Work/Diamond/2021/Al_data/Yun_cameralength_cal/cameralength_cal/Al-standard_0p67_5s/imported.expt",
        refl_file="/home/elena/Work/Diamond/2021/Al_data/Yun_cameralength_cal/cameralength_cal/Al-standard_0p67_5s/strong.refl",
    )

    # params, opt = parse_command_line()
    # Al_data = DialsData(params=params)

    calibrator = PowderCalibrator(Al_data)
    calibrator.tricky_fit_with_calibrant("Al")
