#!/usr/bin/env python
#
# export_json.py
#
#  Copyright (C) 2016 Diamond Light Source
#
#  Author: Richard Gildea
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division

from dials.command_line.reciprocal_lattice_viewer import render_3d

class ReciprocalLatticeJson(render_3d):

  def __init__(self, settings=None):
    render_3d.__init__(self)
    if settings is not None:
      self.settings = settings
    else:
      self.settings = settings()

  def load_models(self, imagesets, reflections):
    self.imagesets = imagesets
    self.reflections_input = reflections
    self.map_points_to_reciprocal_space()

  def as_dict(self, n_digits=None):
    rlp = self.reflections['rlp']
    if n_digits is not None:
      rlp = rlp.round(n_digits)

    if 'imageset_id' in self.reflections:
      imageset_id = list(self.reflections['imageset_id'])
      expt_id = list(self.reflections['id'])
    else:
      imageset_id = list(self.reflections['id'])
      expt_id = None

    indexed = self.reflections.get_flags(self.reflections.flags.indexed)
    d = {
      'rlp': list(rlp),
      'imageset_id': imageset_id,
      'experiment_id': expt_id,
    }
    return d

  def as_json(self, filename=None, compact=False, n_digits=None):
    import json
    d = self.as_dict(n_digits=n_digits)
    if compact:
      text = json.dumps(d, separators=(',',':'), ensure_ascii=True)
    else:
      text = json.dumps(d, indent=2, ensure_ascii=True)

    if filename is not None:
      from libtbx import smart_open
      with smart_open.for_writing(filename) as f:
        f.write(text)
    else:
      return text
