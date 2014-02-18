from __future__ import division

from dials.framework.interface import Interface

class Test(object):

  def __init__(self):
    pass

  def run(self):
    self.tst_before_import()
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
    import dials.interfaces

    # Should have three interfaces
    interfaces = list(Interface.interfaces())
    assert(len(interfaces) == 3)

    # Ensure all the interfaces we expect are there
    from dials.interfaces import CentroidIface
    from dials.interfaces import BackgroundIface
    from dials.interfaces import IntegrationIface

    assert(CentroidIface in interfaces)
    assert(BackgroundIface in interfaces)
    assert(IntegrationIface in interfaces)

    # Should have no extensions
    for iface in interfaces:
      extensions = list(iface.extensions())
      assert(len(extensions) == 0)

    # Test passed
    print 'OK'

  def tst_after_import_extensions(self):
    import dials.extensions
    from dials.interfaces import CentroidIface
    from dials.interfaces import BackgroundIface
    from dials.interfaces import IntegrationIface

    # Should have three interfaces
    interfaces = list(Interface.interfaces())
    assert(len(interfaces) == 3)

    # Check we have the expected number of extensions for each interface
    extensions = list(CentroidIface.extensions())
    assert(len(extensions) == 1)
    extensions = list(BackgroundIface.extensions())
    assert(len(extensions) == 6)
    extensions = list(IntegrationIface.extensions())
    assert(len(extensions) == 4)

    # Check the interface contain the expected extensions
    from dials.extensions import SimpleCentroidExt
    from dials.extensions import NullBackgroundExt
    from dials.extensions import FlatBackgroundExt
    from dials.extensions import InclinedBackgroundExt
    from dials.extensions import CurvedBackgroundExt
    from dials.extensions import FableBackgroundExt
    from dials.extensions import XdsBackgroundExt
    from dials.extensions import Summation2dIntegrationExt
    from dials.extensions import Summation3dIntegrationExt
    from dials.extensions import ProfileFittingRSIntegrationExt
    from dials.extensions import ProfileFittingMosflmIntegrationExt

    extensions = list(CentroidIface.extensions())
    assert(SimpleCentroidExt in extensions)
    extensions = list(BackgroundIface.extensions())
    assert(NullBackgroundExt in extensions)
    assert(FlatBackgroundExt in extensions)
    assert(CurvedBackgroundExt in extensions)
    assert(InclinedBackgroundExt in extensions)
    assert(FableBackgroundExt in extensions)
    assert(XdsBackgroundExt in extensions)
    extensions = list(IntegrationIface.extensions())
    assert(Summation2dIntegrationExt in extensions)
    assert(Summation3dIntegrationExt in extensions)
    assert(ProfileFittingRSIntegrationExt in extensions)
    assert(ProfileFittingMosflmIntegrationExt in extensions)

    # Test passed
    print 'OK'


if __name__ == '__main__':
  test = Test()
  test.run()
