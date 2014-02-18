
class Test(object):

  def __init__(self):
    from dials.framework.config import Config
    self.config = Config()

  def run(self):

    self.config.parse("integration.algorithm=sum3d")

    params = self.config.params()
    assert(params.integration.algorithm == 'sum3d')

    print 'OK'


if __name__ == '__main__':
  test = Test()
  test.run()
