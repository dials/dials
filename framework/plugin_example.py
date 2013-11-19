from __future__ import division
from dials.framework.plugin import Interface, Registry, abstractmethod

class Integrator(object):
  ''' The integration algorithm. '''

  __metaclass__ = Interface

  name = 'integrator'

  @abstractmethod
  def integrate(self):
    pass


class XdsIntegrator(Integrator):
  '''The XDS intensity calculation algorithm.'''

  name = "xds"

  phil = '''
    param1 = 10
      .type = int
      .help = 'a parameter'
  '''

  def __init__(self, sweep, crystal, params):
    super(XdsIntegrator, self).__init__()
    print "Init XDS:"

  def integrate(self):
    print "Integrate"



class MosflmIntegrator(Integrator):
  '''The XDS intensity calculation algorithm.'''

  name = "mosflm"

  def __init__(self, sweep, crystal, params):
    super(MosflmIntegrator, self).__init__()
    print "Init Mosflm:"

  def integrate(self):
    print "Integrate"




class params:
  pass

params.integration = params()
params.integration.algorithm = "xds"
params.integration.xds = params()
params.integration.xds.param1 = 10
params.integration.xds.param2 = 20.5
params.lookup = params()
params.lookup.gain = "A Gain Map"

print list(Registry.extensions(Integrator))

factory = Registry.factory(Integrator)

integrator = factory.create(params.integration.algorithm, None, None, params)
print integrator.integrate()


integrator = factory.create('mosflm', None, None, params)
print integrator.integrate()

print Registry.global_phil().as_str()
