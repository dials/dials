from __future__ import annotations

sacla_phil = """
dispatch.squash_errors = True
dispatch.coset = True
input.reference_geometry=%s
indexing {
  known_symmetry {
    space_group = P43212
    unit_cell = 78.9 78.9 38.1 90 90 90
  }
  refinement_protocol.d_min_start = 2.2
  stills.refine_candidates_with_known_symmetry=True
}
spotfinder {
  filter.min_spot_size = 2
}
refinement {
  parameterisation {
    detector.fix_list = Dist,Tau1
  }
}
profile {
  gaussian_rs {
    centroid_definition = com
  }
}
output.composite_output = True
"""
