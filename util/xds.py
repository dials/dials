from __future__ import annotations

import logging
import os

from dxtbx.serialize import xds
from iotbx.xds import spot_xds
from scitbx import matrix

logger = logging.getLogger(__name__)


def dump(experiments, reflections, directory):
    """Dump the files in XDS format"""
    if len(experiments) > 0:
        for i, experiment in enumerate(experiments):
            suffix = ""
            if len(experiments) > 1:
                suffix = "_%i" % (i + 1)

            sub_dir = f"{directory}{suffix}"
            if not os.path.isdir(sub_dir):
                os.makedirs(sub_dir)
            # XXX imageset is getting the experimental geometry from the image files
            # rather than the input models.expt file
            imageset = experiment.imageset
            imageset.set_detector(experiment.detector)
            imageset.set_beam(experiment.beam)
            imageset.set_goniometer(experiment.goniometer)
            imageset.set_scan(experiment.scan)

            if experiment.crystal is None:
                space_group_number = None
                real_space_a = None
                real_space_b = None
                real_space_c = None
                job_card = "XYCORR INIT COLSPOT IDXREF DEFPIX INTEGRATE CORRECT"
            else:
                crystal_model = experiment.crystal
                crystal_model = crystal_model.change_basis(
                    crystal_model.get_space_group()
                    .info()
                    .change_of_basis_op_to_reference_setting()
                )
                space_group_number = crystal_model.get_space_group().type().number()
                A = matrix.sqr(crystal_model.get_A())
                A_inv = A.inverse()
                real_space_a = A_inv.elems[:3]
                real_space_b = A_inv.elems[3:6]
                real_space_c = A_inv.elems[6:9]
                job_card = ("XYCORR INIT DEFPIX INTEGRATE CORRECT",)

            to_xds = xds.to_xds(imageset)
            xds_inp = os.path.join(sub_dir, "XDS.INP")
            xparm_xds = os.path.join(sub_dir, "XPARM.XDS")
            logger.info("Exporting experiment to %s", xds_inp)
            with open(xds_inp, "w") as f:
                f.write(
                    to_xds.XDS_INP(
                        space_group_number=space_group_number,
                        real_space_a=real_space_a,
                        real_space_b=real_space_b,
                        real_space_c=real_space_c,
                        job_card=job_card,
                    )
                )
            if space_group_number:
                logger.info("Exporting crystal model to %s", xparm_xds)
                with open(xparm_xds, "w") as f:
                    f.write(
                        to_xds.xparm_xds(
                            real_space_a, real_space_b, real_space_c, space_group_number
                        )
                    )

            if reflections is not None and len(reflections) > 0:
                ref_cryst = reflections.select(reflections["id"] == i)
                export_spot_xds(ref_cryst, os.path.join(sub_dir, "SPOT.XDS"))

    else:
        if not os.path.isdir(directory):
            os.makedirs(directory)
        export_spot_xds(reflections, os.path.join(directory, "SPOT.XDS"))


def export_spot_xds(reflections, filename):
    if reflections is not None and len(reflections) > 0:
        centroids = reflections["xyzobs.px.value"]
        intensities = reflections["intensity.sum.value"]
        miller_indices = None
        if "miller_index" in reflections:
            miller_indices = reflections["miller_index"]
            selection = miller_indices != (0, 0, 0)
            miller_indices = miller_indices.select(selection)
            if len(miller_indices) == 0:
                miller_indices = None
            else:
                centroids = centroids.select(selection)
                intensities = intensities.select(selection)
        xds_writer = spot_xds.writer(
            centroids=centroids, intensities=intensities, miller_indices=miller_indices
        )
        logger.info("Exporting spot list as %s", filename)
        xds_writer.write_file(filename=filename)
