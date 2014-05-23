from __future__ import division
# script to convert the output of refined_quadrants to a cbf header that
# can be applied to a cspad cbf with cxi.apply_metrology
# Note the hardcoded wavelength and distance.  For development only.

from dials.util.command_line import Importer

importer = Importer(['refined_experiments.json'], check_format=False)
experiment = importer.experiments[0]
detector = experiment.detector

from xfel.cftbx.detector.cspad_cbf_tbx import write_cspad_cbf, map_detector_to_basis_dict

metro = map_detector_to_basis_dict(detector)
write_cspad_cbf(None, metro, 'cbf', None, 'quad_refined.def', 1.81587, 105, header_only=True)

print "Done"
