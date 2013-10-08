
from dials.framework.toplevel import Interface

class HighLevelIntegrator(Interface):

    def __init__(self):
        super(HighLevelIntegrator, self).__init__()

    def prepare(self):
        self.prepared = True
        print "Prepare:"
        print " Predict reflections"
        print " Extract shoeboxes"

    def process(self):
        self.processed = True
        print "Do:"
        print " Calculate centroids"
        print " Calculate background"
        print " Calculate intensities"

    def finish(self):
        self.finished = True
        print "Finished:"
        print " Perform corrections"


integrator = HighLevelIntegrator()
integrator.run()
