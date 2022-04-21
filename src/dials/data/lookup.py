from __future__ import annotations

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
