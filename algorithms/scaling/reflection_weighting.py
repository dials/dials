from dials.array_family import flex

class Weighting(object):
  ''' This class defines a weighting object that takes in a reflection table,
  gives initial weights of 1/variance and has methods to set the weights of
  certain reflections to zero.'''
  def __init__(self, reflection_table):
    '''set initial weighting to be a statistical weighting'''
    self.scale_weighting = 1.0/reflection_table['variance']

  def get_weights(self):
    '''method to return flex array containing weights'''
    return self.scale_weighting

  def set_unity_weighting(self, reflection_table):
    '''method to weight each reflection equally'''
    self.scale_weighting = flex.double([1.0]*len(reflection_table['variance']))

  def apply_Isigma_cutoff(self, reflection_table, ratio):
    '''method to set a zero weight below an I/sigma cutoff'''
    Ioversigma = reflection_table['intensity']/(reflection_table['variance']**0.5)
    sel = Ioversigma <= ratio
    self.scale_weighting.set_selected(sel, 0.0)

  def apply_dmin_cutoff(self, reflection_table, d_cutoff_value):
    '''method to set a zero weight below an d-value cutoff'''
    sel = reflection_table['d'] <= d_cutoff_value
    self.scale_weighting.set_selected(sel, 0.0)

  def remove_wilson_outliers(self, reflection_table):
    '''method to set a zero weight for any outliers as determined by
    Wilson_outlier_test.py'''
    if 'wilson_outlier_flag' in reflection_table:
      sel = reflection_table['wilson_outlier_flag']
      self.scale_weighting.set_selected(sel, 0.0)

  def apply_aimless_error_model(self, reflection_table, error_params):
    '''def new_sigma(weight):
      sigmaprime = error_params[0] * (((1.0/weight)
                                     + ((error_params[1] * reflection_table['intensity'])**2))**0.5)
      return sigmaprime
    new_weights = [1.0/(new_sigma(i)**2) if i > 0.0 else 0.0 for i in self.scale_weighting]'''
    sigmaprime = error_params[0] * (((1.0/self.scale_weighting)
                                     + ((error_params[1] * reflection_table['intensity'])**2))**0.5)
    self.scale_weighting = 1.0/(sigmaprime**2)
