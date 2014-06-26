from __future__ import division

class Test(object):

  def __init__(self):
    from dials.framework.registry import Registry
    self.config = Registry().config()

  def run(self):

    self.config.parse("integration.intensity.algorithm=sum3d")

    params = self.config.params()
    assert(params.integration.intensity.algorithm == 'sum3d')

    print 'OK'


if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    test = Test()
    test.run()
