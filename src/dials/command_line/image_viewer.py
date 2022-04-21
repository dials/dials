# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1


from __future__ import annotations

import pickle
import sys

import iotbx.phil

import dials.util.log
from dials.util.image_viewer.spotfinder_wrap import spot_wrapper
from dials.util.options import ArgumentParser, flatten_experiments, flatten_reflections

help_message = """

This program can be used for viewing diffraction images, optionally overlaid
with the results of spot finding, indexing or integration.

Examples::

  dials.image_viewer image.cbf

  dials.image_viewer models.expt

  dials.image_viewer models.expt strong.refl

  dials.image_viewer models.expt integrated.refl
"""

phil_scope = iotbx.phil.parse(
    """\
brightness = 100
  .type = int
color_scheme = *grayscale rainbow heatmap invert
  .type = choice
projection = lab *image
  .type = choice
show_beam_center = True
  .type = bool
show_resolution_rings = False
  .type = bool
show_ice_rings = False
  .type = bool
show_ctr_mass = True
  .type = bool
show_max_pix = True
  .type = bool
show_all_pix = True
  .type = bool
show_threshold_pix = False
  .type = bool
show_shoebox = True
  .type = bool
show_predictions = True
  .type = bool
show_miller_indices = False
  .type = bool
show_indexed = False
  .type = bool
show_integrated = False
  .type = bool
show_mask = False
  .type = bool
show_basis_vectors = True
  .type = bool
show_rotation_axis = False
  .type = bool
basis_vector_scale = 10
  .type = int(value_min=1, value_max=20)
display = *image mean variance dispersion sigma_b \
          sigma_s threshold global_threshold
  .type = choice
nsigma_b = 6
  .type = float(value_min=0)
nsigma_s = 3
  .type = float(value_min=0)
global_threshold = 0
  .type = float(value_min=0)
kernel_size = 3,3
  .type = ints(size=2, value_min=1)
min_local = 2
  .type = int
gain = 1
  .type = float(value_min=0)
  .help = "Set gain for the dispersion algorithm. This does not override the"
          "detector's panel gain, but acts as a multiplier for it."

include scope dials.extensions.radial_profile_spotfinder_threshold_ext.phil_str

stack_images = 1
  .type = int(value_min=1)
  .expert_level = 2
stack_mode = max mean *sum
  .type = choice
d_min = None
  .type = float(value_min=0)
mask = None
  .type = path
  .help = path to mask pickle file

powder_arcs{
  show = False
    .type=bool
    .help = "show powder arcs calculated from PDB file."
  code = None
    .type=str
    .help = "PDB code (4 characters) for file; fetch it from the Internet."
}
calibrate_silver = False
    .type=bool
    .help = "Open special GUI for distance/metrology from silver behenate."
calibrate_pdb{
  code = None
    .type=str
    .help = "Open pdb code (over Internet) to get unit cell & symmetry for powder rings."
    .help = "Most useful for calibrating low-Q rings on far detector."
    .help = "Option is mutually exclusive with calibrate silver, unit cell and powder arcs options."
  d_min = 20.
    .type=float
    .help = "Limiting resolution to calculate powder rings"
}
calibrate_unit_cell{
  unit_cell = None
    .type=unit_cell
    .help = "Specify unit cell for powder rings."
    .help = "Option is mutually exclusive with calibrate silver, pdb and powder arcs options."
  d_min = 20.
    .type=float
    .help = "Limiting resolution to calculate powder rings"
  space_group = None
    .type=str
    .help = "Specify spacegroup for the unit cell"
  show_hkl = None
    .type = ints(size=3)
    .multiple=True
    .help = "Limit display of rings to these Miller indices"
}

include scope dials.util.options.format_phil_scope

masking {
  include scope dials.util.masking.phil_scope
}

output {
  mask = pixels.mask
    .type = path
    .help = "Name of output mask file"

  mask_params = mask.phil
    .type = path
    .help = "Name of output mask parameter file"
}

predict_reflections = False
  .type = bool
  .help = Predict reflections if no reflections provided in input

include scope dials.algorithms.profile_model.factory.phil_scope
include scope dials.algorithms.spot_prediction.reflection_predictor.phil_scope

load_models = True
  .type = bool
  .help = "Whether to load every model, which matters for large image files"

zmq_endpoint = None
  .type = str
  .help = "The endpoint to bind a zeromq PULL socket to, for receiving commands"
  .expert_level = 3
""",
    process_includes=True,
)


def show_image_viewer(params, experiments, reflections):
    wrapper = spot_wrapper(params=params)
    wrapper.display(experiments=experiments, reflections=reflections)


@dials.util.show_mail_handle_errors()
def run(args=None):
    dials.util.log.print_banner()
    usage_message = "dials.image_viewer models.expt [observations.refl]"
    parser = ArgumentParser(
        usage=usage_message,
        phil=phil_scope,
        read_experiments=True,
        read_reflections=True,
        read_experiments_from_images=True,
        epilog=help_message,
    )
    params, options = parser.parse_args(args, show_diff_phil=True)
    experiments = [x.data for x in params.input.experiments]
    reflections = flatten_reflections(params.input.reflections)

    if len(experiments) == 0:
        parser.print_help()
        exit(0)

    flat_expts = flatten_experiments(params.input.experiments)
    if params.load_models:
        if any(e.detector is None for e in flat_expts):
            sys.exit("Error: experiment has no detector")
        if any(e.beam is None for e in flat_expts):
            sys.exit("Error: experiment has no beam")

    # If given a mask, replace the path with the loaded data
    if params.mask is not None:
        with open(params.mask, "rb") as f:
            params.mask = pickle.load(f)

    show_image_viewer(params=params, reflections=reflections, experiments=experiments)


if __name__ == "__main__":
    run()
