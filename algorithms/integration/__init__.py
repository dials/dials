from __future__ import division
import boost.python
from cctbx import sgtbx # import dependency
from dials.array_family import flex # import dependency
from dials.model.data import AdjacencyList # import dependency
from dials_algorithms_integration_ext import *
from integrator import *
from integrator_stills import *


class ImageSummaryAux(boost.python.injector, ImageSummary):

  def table(self):
    from libtbx.table_utils import format as table
    rows = [["Image",
             "# full",
             "# part",
             "<I/sigI>\n (sum)",
             "<I/sigI>\n (prf)"]]
    for (i,
         (full,
          part,
          sum_ios,
          prf_ios)) in enumerate(zip(
           self.full(),
           self.part(),
           self.sum_ios(),
           self.prf_ios())):
      rows.append([
        '%d' % i,
        '%d' % full,
        '%d' % part,
        '%.1f' % sum_ios,
        '%.1f' % prf_ios])
    return table(rows, has_header=True, justify='right', prefix=' ')

class ResolutionSummaryAux(boost.python.injector, ResolutionSummary):

  def table(self):
    from libtbx.table_utils import format as table
    d = self.bins()
    rows = [["d",
             "<I/sigI>\n (sum)",
             "<I/sigI>\n (prf)"]]
    for (i,
         (sum_ios,
          prf_ios)) in enumerate(zip(
           self.sum_ios(),
           self.prf_ios())):
      rows.append([
        '%.1f -> %.1f' % (d[i], d[i+1]),
        '%.1f' % sum_ios,
        '%.1f' % prf_ios])
    return table(rows, has_header=True, justify='right', prefix=' ')

class WholeSummaryAux(boost.python.injector, WholeSummary):

  def table(self):
    from libtbx.table_utils import format as table
    rows = [["<I/sigI>\n (sum)",
             "<I/sigI>\n (prf)"]]
    rows.append([
      '%.1f' % self.sum_ios(),
      '%.1f' % self.prf_ios()])
    return table(rows, has_header=True, justify='right', prefix=' ')

class SummaryAux(boost.python.injector, Summary):

  def as_str(self):
    img_summary = self.image_summary().table()
    res_summary = self.resolution_summary().table()
    who_summary = self.whole_summary().table()
    return (
      ' Summary of integration results as a function of image number'
      '\n%s\n\n'
      ' Summary of integration results binned by resolution'
      '\n%s\n\n'
      ' Summary of integration results for the whole dataset'
      '\n%s\n'
    ) % (img_summary, res_summary, who_summary)
