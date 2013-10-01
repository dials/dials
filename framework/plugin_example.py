from dials.framework import plugin

class Integrator(object):
    '''The intensity calculation algorithm.'''
    __metaclass__ = plugin.Extension

    name = "integration"


class XdsPlugin(Integrator):
    '''The XDS intensity calculation algorithm.'''

    name = "xds"

    config = '''
        param1 = 10
          .type = int
          .help = "This is a parameter"
        param2 = 20.5
          .type = float
          .help = "This is another parameter"
    '''

    def __init__(self, sweep, crystal, params):
        print "Init XDS:"
        print "  ", sweep, crystal
        print "  ", params.integration.xds.param1
        print "  ", params.integration.xds.param2
        print "  ", params.lookup.gain

    def __call__(self):
        print "Call XDS"


class MosflmPlugin(Integrator):
    '''The XDS intensity calculation algorithm.'''

    name = "mosflm"

    def __init__(self, sweep, crystal, params):
        print "Init Mosflm"

    def __call__(self):
        print "Call Mosflm"



extensions = plugin.extensions()
print extensions
print extensions.plugins()
print extensions.configuration()
print Integrator.plugins()
print Integrator.configuration()

class params:
    pass

params.integration = params()
params.integration.algorithm = "xds"
params.integration.xds = params()
params.integration.xds.param1 = 10
params.integration.xds.param2 = 20.5
params.lookup = params()
params.lookup.gain = "A Gain Map"


calculate_intensity = Integrator.factory(params.integration.algorithm, None, None, params)
calculate_intensity()

calculate_intensity = Integrator.factory("mosflm", None, None, params)
calculate_intensity()
