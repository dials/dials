from __future__ import division
from cctbx.array_family import flex # import dependency
from dials_model_data_ext import *

def reflection_list_getinitargs(self):
    ''' Pickle the reflection list by converting to a list. '''
    return (list(self),)

# Inject the __getinitargs__ method into ReflectionList
ReflectionList.__getinitargs__ = reflection_list_getinitargs
