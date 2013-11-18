from __future__ import division
from scitbx import matrix

class source:
  """Simple model for the source X-ray beam in which the beam vector can be
  moved, but changes to its initial magnitude are not allowed"""

  def __init__(self, s0_vec, wavelength = None):
    """Beam direction is taken from s0_vec. If wavelength is supplied that
    sets the magnitude, otherwise the length of s0_vec is preserved"""
    assert isinstance(s0_vec, matrix.rec)
    self._s0 = matrix.col(s0_vec)
    if wavelength is not None:
      self._s0 = self._s0.normalize() / wavelength

  def get_s0(self):
    """Return the reciprocal space beam vector"""
    return self._s0

  def get_beam(self):
    """Return the real space beam vector"""
    return self._s0 * 1. / self._s0.length_sq()

  def get_wavelength(self):
    """Return the wavelength"""
    return 1. / self._s0.length()

  def set_s0(self, vals):
    """set the beam vector, but don't allow changes to the wavelength"""
    if len(vals) != 3:
      assert isinstance(vals, matrix.col)
    new_s0 = matrix.col(vals)
    mag_ratio = self._s0.length() / new_s0.length()
    self._s0 = mag_ratio * new_s0
