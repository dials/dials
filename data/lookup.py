from __future__ import absolute_import, division, print_function

from libtbx.phil import parse

phil_scope = parse(
    """

lookup
  .help = "Parameters specifying lookup file path"
{
  mask = None
    .help = "The path to the mask file."
    .type = str
}
"""
)
