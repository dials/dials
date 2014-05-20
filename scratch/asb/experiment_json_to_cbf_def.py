from __future__ import division

# Script to convert the output from refine_quadrants to a header file
# Apply the header file to the cbfs with cxi.apply_metrology
# Note hardcoded distance of 105

from dials.util.command_line import Importer
from xfel.cftbx.detector.cspad_cbf_tbx import write_cspad_cbf, map_detector_to_basis_dict

importer = Importer(['refined_experiments.json'], check_format=False)
experiment = importer.experiments[0]
detector = experiment.detector

metro = map_detector_to_basis_dict(detector)
write_cspad_cbf(None, metro, 'cbf', None, 'quad_refined.def', None, 105, header_only=True)

print "Done"
