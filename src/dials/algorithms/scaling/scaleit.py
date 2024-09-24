mtz1 = "merged_3.mtz"
mtz2 = "merged_1.mtz"
expts = "scaled.expt"

anomalous = False

from dials.array_family import flex
from dxtbx.model import Experiment, ExperimentList
from dials.util.reference import intensities_from_reference_data_file
from iotbx.mtz.crystal_symmetry_from_mtz import extract_from
from dxtbx.model import Crystal
from iotbx import mtz
from dxtbx.serialize import load
elist = load.experiment_list(expts)
from copy import deepcopy
from unittest.mock import Mock
from dials.algorithms.scaling.model.model import KBScalingModel
from dxtbx.util import ersatz_uuid4

from dials.util.options import ArgumentParser
from libtbx import Auto, phil
import logging
from dials.util import log
log.config()
logger = logging.getLogger("dials")

def datafiles_from_mtz(mtzfile, id_=0):
    intensities = intensities_from_reference_data_file(mtzfile)
    if anomalous:
        intensities = intensities.as_non_anomalous_array().merge_equivalents().array()
    table = flex.reflection_table()
    table["intensity"] = intensities.data()
    table["miller_index"] = intensities.indices()
    table["id"] = flex.int(table.size(), id_)
    table["d"] = intensities.d_spacings().data()
    table["variance"] = intensities.sigmas() ** 2
    table.set_flags(flex.bool(table.size(), False), table.flags.bad_for_scaling)
    table.set_flags(flex.bool(table.size(), True), table.flags.integrated)

    expt = Experiment()
    expt.crystal = deepcopy(elist[0].crystal)
    phil_scope = phil.parse(
        """
    include scope dials.algorithms.scaling.scaling_options.phil_scope
    include scope dials.algorithms.scaling.model.model.model_phil_scope
    include scope dials.algorithms.scaling.scaling_refiner.scaling_refinery_phil_scope
""",
        process_includes=True,
    )
    parser = ArgumentParser(phil=phil_scope, check_format=False)
    params, _ = parser.parse_args(args=[], quick_parse=True)
    params.model = "KB"
    params.KB.decay_correction = False
    expt.scaling_model = KBScalingModel.from_data(params, [], [])
    #expt.scaling_model.set_scaling_model_as_scaled()  # Set as scaled to fix scale.
    expt.identifier = ersatz_uuid4()
    table.experiment_identifiers()[id_] = expt.identifier
    return table, expt

t1,e1 = datafiles_from_mtz(mtz1, 0)
e1 = ExperimentList([e1])
t2,e2 = datafiles_from_mtz(mtz2, 1)
e2 = ExperimentList([e2])
t1["intensity.sum.value"] = t1["intensity"]
t1["intensity.sum.variance"] = t1["variance"]
t2["intensity.sum.value"] = t2["intensity"]
t2["intensity.sum.variance"] = t2["variance"]
from dials.algorithms.scaling.scaling_library import scale_against_target

result = scale_against_target(t1,e1,t2,e2)
print(e1.scaling_models()[0].to_dict())
print(e1.scaling_models()[1].to_dict())


print("done")