class MosflmProfileFitting:
  def __init__(self, nblocks):

    self.nblocks = nblocks
  def __call__(self, sweep, crystal, rlist, reference = None):
    old_way='''
    from dials.algorithms.integration import flex_2d_layering_n_integrating
    from dials.scratch.luiso_s.test_code.call_mosflm_2d import mosflm_caller

    from dials.algorithms.background.curved_background_subtractor \
     import  layering_and_background_modl

    layering_and_background_modl(rlist)
    flex_2d_layering_n_integrating(rlist)
    xmax, ymax = sweep.get_detector()[0].get_image_size()
    rlist = mosflm_caller(rlist, xmax, ymax, self.nblocks)
    return rlist
    '''

    coorected_old_way='''
    from dials.algorithms.integration import flex_2d_layering_n_integrating
    from dials.scratch.luiso_s.test_code.call_mosflm_2d import mosflm_caller

    #from dials.algorithms.background.curved_background_subtractor \
    # import  layering_and_background_modl

    #layering_and_background_modl(rlist)
    flex_2d_layering_n_integrating(rlist)
    xmax, ymax = sweep.get_detector()[0].get_image_size()
    rlist = mosflm_caller(rlist, xmax, ymax, self.nblocks)
    '''

    from dials.algorithms.integration import flex_2d_layering_n_integrating
    #from dials.scratch.luiso_s.test_code.call_mosflm_2d import mosflm_caller
    from dials.algorithms.integration.call_mosflm_2d  import mosflm_caller

    flex_2d_layering_n_integrating(rlist)
    xmax, ymax = sweep.get_detector()[0].get_image_size()
    rlist = mosflm_caller(rlist, xmax, ymax, self.nblocks)
