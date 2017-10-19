from dials.array_family import flex
import numpy as np

'''
This class defines a weighting object that takes in a reflection table, gives
initial weights of 1/variance and has methods to set the weights of certain
reflections to zero.
'''

class Weighting(object):
  def __init__(self, reflection_table):
    '''set initial weighting to be a statistical weighting'''
    self.scale_weighting = 1.0/reflection_table['variance']

  def get_weights(self):
    return self.scale_weighting

  def set_unity_weighting(self, reflection_table):
    self.scale_weighting = flex.double([1.0]*len(reflection_table['variance']))

  def apply_Isigma_cutoff(self, reflection_table, ratio):
    sel = flex.bool()
    for i, intensity in enumerate(reflection_table['intensity']):
      if ratio > intensity/(reflection_table['variance'][i]**0.5):
        sel.append(True)
      else:
        sel.append(False)
    self.scale_weighting.set_selected(sel, 0.0)
    #print len(self.scale_weighting) - self.scale_weighting.count(0.0)

  def apply_dmin_cutoff(self, reflection_table, d_cutoff_value):
    sel = flex.bool()
    for i, d in enumerate(reflection_table['d']):
      if d_cutoff_value > d:
        sel.append(True)
      else:
        sel.append(False)
    self.scale_weighting.set_selected(sel, 0.0)
    #print len(self.scale_weighting) - self.scale_weighting.count(0.0)
  
  def remove_wilson_outliers(self, reflection_table):
    if 'wilson_outlier_flag' in reflection_table:
      sel = reflection_table['wilson_outlier_flag']
      self.scale_weighting.set_selected(sel, 0.0)