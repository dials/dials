from __future__ import division

from dials.framework.interface import Interface

class Test(object):

  def __init__(self):
    pass

  def run(self):
    # self.tst_before_import()
    self.tst_after_import_interfaces()
    self.tst_after_import_extensions()

  def tst_before_import(self):

    # Currently no interfaces
    interfaces = list(Interface.interfaces())
    assert(len(interfaces) == 0)

    # Interface class has no extensions
    try:
      Inteface.extensions()
      assert(False)
    except Exception:
      pass

    # Test passed
    print 'OK'

  def tst_after_import_interfaces(self):
    import dials.interfaces # import dependency

    # Should have three interfaces
    interfaces = list(Interface.interfaces())
    assert(len(interfaces) == 4)

    # Ensure all the interfaces we expect are there
    from dials.interfaces import ProfileModelIface
    from dials.interfaces import SpotFinderThresholdIface
    from dials.interfaces import CentroidIface
    from dials.interfaces import BackgroundIface

    assert(ProfileModelIface in interfaces)
    assert(SpotFinderThresholdIface in interfaces)
    assert(CentroidIface in interfaces)
    assert(BackgroundIface in interfaces)

    # Test passed
    print 'OK'

  def tst_after_import_extensions(self):
    import dials.extensions # import dependency
    from dials.interfaces import ProfileModelIface
    from dials.interfaces import SpotFinderThresholdIface
    from dials.interfaces import CentroidIface
    from dials.interfaces import BackgroundIface

    # Should have four interfaces
    interfaces = list(Interface.interfaces())
    assert(len(interfaces) == 4)

    # Check we have the expected number of extensions for each interface
    extensions = list(ProfileModelIface.extensions())
    assert(len(extensions) > 0)
    extensions = list(SpotFinderThresholdIface.extensions())
    assert(len(extensions) > 0)
    extensions = list(CentroidIface.extensions())
    assert(len(extensions) > 0)
    extensions = list(BackgroundIface.extensions())
    assert(len(extensions) > 0)

    # Check the interface contain the expected extensions
    from dials.extensions import GaussianRSProfileModelExt
    from dials.extensions import KabschSpotFinderThresholdExt
    from dials.extensions import SimpleCentroidExt
    from dials.extensions import NullBackgroundExt
    from dials.extensions import SimpleBackgroundExt

    extensions = list(ProfileModelIface.extensions())
    assert(GaussianRSProfileModelExt in extensions)
    extensions = list(SpotFinderThresholdIface.extensions())
    assert(KabschSpotFinderThresholdExt in extensions)
    extensions = list(CentroidIface.extensions())
    assert(SimpleCentroidExt in extensions)
    extensions = list(BackgroundIface.extensions())
    assert(NullBackgroundExt in extensions)
    assert(SimpleBackgroundExt in extensions)

    # Test passed
    print 'OK'


if __name__ == '__main__':
  from dials.test import cd_auto
  try:
    with cd_auto(__file__):
      test = Test()
      test.run()
  except Exception, e:

    message = '''
    ====================================================
      Try deleting old *.pyc files from dials/extensions
    ====================================================
    '''

    e.args = (str(e.args) + '\n' + message,)
    raise
