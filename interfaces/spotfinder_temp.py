from dials.framework import toplevel

class SpotFinderInterface(toplevel.Interface):

    def __init__(self, sweep, crystal, params):

        from dials.framework.plugin import Registry
        from dials.interfaces.standard_interfaces import SpotFinderThreshold
        from dials.interfaces.standard_interfaces import SpotFinderFilter

        threshold_factory = Registry.factory(SpotFinderThreshold)
        filter_factory = Registry.factory(SpotFinderFilter)

    def prepare(self):
        self.prepared = True

    def process(self):

        self.processed = True

    def finish(self):
        self.finished = True
