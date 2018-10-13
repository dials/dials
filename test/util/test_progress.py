"""
Tests for the dials-customized version of the tqdm progress bar
"""

from __future__ import absolute_import, division, print_function

import time as systime

import pytest
from six.moves import StringIO

from dials.util import progress
import tqdm


class TTYStringIO(StringIO):
  """A StringIO Subclass that states it is a TTY"""

  def isatty(self):
    return True

  def getvalue(self):
    # Remove escape characters so that the output makes sense if printed
    return StringIO.getvalue(self).replace("\x1b", r"\x1b")

  def write(self, x):
    return StringIO.write(self, x)


class FakeTime(object):
  """Fake time class that lets us control what has passed"""

  def __init__(self):
    self._time = systime.time()

  def __call__(self):
    return self._time

  def sleep(self, elapsed):
    self._time += elapsed

  def time(self):
    return self._time


@pytest.fixture
def bar(time):
  """Simple, repeating setup for a progress bar for tests"""
  io = TTYStringIO()
  pbar = progress(file=io, total=10)
  pbar.io = io
  # Validate the time module was intercepted
  assert pbar._time is time
  yield pbar
  pbar.close()


@pytest.fixture
def time(monkeypatch):
  "Fixture to monkeypatch the tqdm time access"
  ft = FakeTime()
  monkeypatch.setattr(tqdm._tqdm, "time", ft)
  return ft


def test_noshow(bar, time):
  "Basic test of not showing the bar"
  # Make sure that if the ETA stays under 10 we never show anything
  for x in range(10):
    time.sleep(0.5)
    bar.update(1)
    assert not bar.io.getvalue()


def test_show(bar, time):
  """Test that the bar shows for an ETA over 10s"""
  time.sleep(2.2)
  bar.update(2)
  assert bar.io.getvalue()


def test_notty(time):
  "Test that nothing shows for a non-tty"
  # Setup a bar with a non-tty file
  io = StringIO()
  bar = progress(file=io, total=10)
  # Take steps that would normally show
  time.sleep(11)
  bar.update(1)
  time.sleep(100)
  bar.update(9)
  assert not io.getvalue()
  bar.close()


def test_slowstart(bar, time):
  "Test the case where a slow first entry lies"
  assert not bar.io.getvalue()
  bar.total = 100
  # Give an update that estimates 20 seconds
  time.sleep(0.2)
  bar.update(1)
  assert not bar.io.getvalue()
  for x in range(99):
    time.sleep(0.05)
    bar.update(1)
    assert not bar.io.getvalue()
