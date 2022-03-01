from __future__ import annotations

from libtbx.phil import parse

phil_scope = parse(
    """

mp {
  method = *multiprocessing sge lsf pbs
    .type = choice
    .help = "The multiprocessing method to use"

  nproc = 1
    .type = int(value_min=1)
    .help = "The number of processes to use."

  nthreads = 1
    .type = int(value_min=1)
    .help = "The number of local threads to use for openmp."
}
"""
)
