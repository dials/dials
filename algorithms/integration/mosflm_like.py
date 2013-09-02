class MosflmProfileFitting:
    def __init__(self, nblocks):

        self.nblocks = nblocks
    def __call__(self, sweep, crystal, rlist, reference = None):
        from dials.algorithms.integration import flex_2d_layering_n_integrating
        from dials.scratch.luiso_s.test_code.call_mosflm_2d import mosflm_caller

        flex_2d_layering_n_integrating(rlist)
        xmax, ymax = sweep.get_detector().get_image_size()
        rlist = mosflm_caller(rlist, xmax, ymax, self.nblocks)
        return rlist

