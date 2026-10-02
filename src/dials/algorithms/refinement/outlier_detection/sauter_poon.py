from __future__ import annotations

from dials.algorithms.refinement.outlier_detection import CentroidOutlier


class SauterPoon(CentroidOutlier):
    """Implementation of the CentroidOutlier class using the algorithm of
    Sauter & Poon (2010) (https://doi.org/10.1107/S0021889810010782)."""

    def __init__(
        self,
        cols=None,
        min_num_obs=20,
        separate_experiments=True,
        separate_panels=True,
        separate_images=False,
        block_width=None,
        nproc=1,
        px_sz=(1, 1),
        verbose=False,
        pdf=None,
    ):
        # here the column names are fixed by the algorithm, so what's passed in is
        # ignored.
        CentroidOutlier.__init__(
            self,
            cols=["miller_index", "xyzobs.px.value", "xyzcal.px", "panel"],
            min_num_obs=min_num_obs,
            separate_experiments=separate_experiments,
            separate_panels=separate_panels,
            separate_images=separate_images,
            block_width=block_width,
            nproc=nproc,
        )

        # px_sz is either a single (x, y) pixel size in mm applied to all panels,
        # or a sequence of (x, y) pixel sizes, one per panel, for detectors
        # whose panels have different pixel sizes
        self._px_sz = px_sz
        self._verbose = verbose
        self._pdf = pdf

        return

    def _px_sz_for_panel(self, ipanel):
        try:
            return self._px_sz[ipanel]
        except (TypeError, IndexError):
            return self._px_sz

    def _detect_outliers(self, cols):
        # cols is guaranteed to be a list of four flex arrays, containing miller
        # indices, observed pixel coordinates, calculated pixel coordinates and
        # panel ids. Copy the data into matches
        class match:
            pass

        single_px_sz = len(self._px_sz) == 2 and not hasattr(self._px_sz[0], "__len__")

        matches = []
        for hkl, obs, calc, ipanel in zip(cols[0], cols[1], cols[2], cols[3]):
            m = match()
            m.miller_index = hkl
            px_sz = self._px_sz if single_px_sz else self._px_sz_for_panel(ipanel)
            m.x_obs = obs[0] * px_sz[0]
            m.y_obs = obs[1] * px_sz[1]
            m.x_calc = calc[0] * px_sz[0]
            m.y_calc = calc[1] * px_sz[1]
            matches.append(m)

        import iotbx.phil
        from rstbx.phil.phil_preferences import indexing_api_defs

        hardcoded_phil = iotbx.phil.parse(input_string=indexing_api_defs).extract()

        # set params into the hardcoded_phil
        hardcoded_phil.indexing.outlier_detection.verbose = self._verbose
        hardcoded_phil.indexing.outlier_detection.pdf = self._pdf

        from rstbx.indexing_api.outlier_procedure import OutlierPlotPDF

        if self._pdf is not None:
            ## new code for outlier rejection inline here
            hardcoded_phil.__inject__(
                "writer", OutlierPlotPDF(hardcoded_phil.indexing.outlier_detection.pdf)
            )

        # execute Sauter and Poon (2010) algorithm
        from rstbx.indexing_api import outlier_detection

        od = outlier_detection.find_outliers_from_matches(
            matches, verbose=self._verbose, horizon_phil=hardcoded_phil
        )

        # flex.bool of the inliers
        outliers = ~od.get_cache_status()

        return outliers
