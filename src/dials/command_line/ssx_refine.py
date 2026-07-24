# DIALS_ENABLE_COMMAND_LINE_COMPLETION
"""
This program runs the dials.refine program with suitable defaults for
joint detector geometry refinement of ssx data. See the dials.refine
documentation for full details of the program.

Examples::

  dials.ssx_refine indexed.expt indexed.refl

"""

from __future__ import annotations

import libtbx.phil

from dials.command_line.refine import run as refine_run
from dials.command_line.refine import working_phil as dials_refine_phil

usage = "dials.ssx_refine indexed.expt indexed.refl [options] [param.phil] "

# local overrides for refiner.phil_scope
phil_overrides = libtbx.phil.parse(
    """
refinement {
  parameterisation {
    auto_reduction {
      action = fail *fix remove
    }
    beam {
      fix = *all in_spindle_plane out_spindle_plane wavelength
    }
    detector {
      fix_list = "Tau1"
    }
  }
  refinery {
    engine = SimpleLBFGS LBFGScurvs GaussNewton LevMar *SparseLevMar
  }
  reflections {
    outlier {
      algorithm = null auto mcd tukey *sauter_poon
    }
  }
}
output {
    log = "dials.ssx_refine.log"
}
"""
)

working_phil = dials_refine_phil.fetch(sources=[phil_overrides])


def run(args=None):
    refine_run(args=args, phil=working_phil, usage=usage, epilog=__doc__)


if __name__ == "__main__":
    run()
