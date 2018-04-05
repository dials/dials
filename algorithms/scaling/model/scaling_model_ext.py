#!/usr/bin/env python
"""
Definitions of scaling model extensions.
"""
from __future__ import absolute_import, division

from dials.framework import interface

class ScalingModelIface(interface.Interface):
    '''
    The interface definition for a scaling model.
    '''
    scope = "scaling"
    name = 'scaling_model'

    @classmethod
    def factory(cls):
        ''' Get the factory. '''
        pass

    @classmethod
    def scaler(cls):
        ''' Get the scaler. '''
        pass

    @staticmethod
    def from_dict(d):
        ''' Get from dictionary. '''
        pass

class PhysicalScalingModelExt(ScalingModelIface):
  """An extension class implementing a physical scaling model."""

  name = 'physical'

  @classmethod
  def factory(cls):
    '''returns the scaling Model Factory'''
    from dials.algorithms.scaling.model.scaling_model_factory import \
      PhysicalSMFactory
    return PhysicalSMFactory

  @classmethod
  def scaler(cls):
    '''returns the scaler factory'''
    from dials.algorithms.scaling.scaler import PhysicalScaler
    return PhysicalScaler

  @staticmethod
  def from_dict(d):
    '''creates a scaling model from a dict'''
    from dials.algorithms.scaling.model.model import \
      PhysicalScalingModel
    return PhysicalScalingModel.from_dict(d)

class KBScalingModelExt(ScalingModelIface):
  """An extension class implementing a KB scaling model."""

  name = 'KB'

  @classmethod
  def factory(cls):
    '''returns the scaling Model Factory'''
    from dials.algorithms.scaling.model.scaling_model_factory import \
      KBSMFactory
    return KBSMFactory

  @classmethod
  def scaler(cls):
    '''returns the scaler factory'''
    from dials.algorithms.scaling.scaler import KBScaler
    return KBScaler

  @staticmethod
  def from_dict(d):
    '''creates a scaling model from a dict'''
    from dials.algorithms.scaling.model.model import KBScalingModel
    return KBScalingModel.from_dict(d)

class ArrayScalingModelExt(ScalingModelIface):
  """An extension class implementing an array-based scaling model."""

  name = 'array'

  @classmethod
  def factory(cls):
    '''returns the scaling Model Factory'''
    from dials.algorithms.scaling.model.scaling_model_factory import \
      ArraySMFactory
    return ArraySMFactory

  @classmethod
  def scaler(cls):
    '''returns the scaler factory'''
    from dials.algorithms.scaling.scaler import ArrayScaler
    return ArrayScaler

  @staticmethod
  def from_dict(d):
    '''creates a scaling model from a dict'''
    from dials.algorithms.scaling.model.model import ArrayScalingModel
    return ArrayScalingModel.from_dict(d)
