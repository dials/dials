from __future__ import absolute_import, division

class Test(object):
  def run(self):
    from libtbx import easy_run

    # Call dials.merge_reflection_lists
    easy_run.fully_buffered([
      'dev.dials.show_extensions',
    ]).raise_if_errors()

    print 'OK'

if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    test = Test()
    test.run()
