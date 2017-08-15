import pickle as pickle
from dials_array_family_flex_ext import *
from cctbx.array_family.flex import *
from cctbx.array_family import flex
from cctbx import miller
from cctbx import crystal

from dials.util.options import Importer, flatten_experiments, flatten_reflections, OptionParser




class Data_Manager:
    def __init__(self,reflection_table,experimental_unit_cell,experimental_space_group):
        self.reflection_table = reflection_table
        xyzvalues = reflection_table['xyzobs.px.value']
        zvalues = flex.double([xyzvalues[i][2] for i in range(0,len(xyzvalues))])
        self.reflection_table['z_value'] = zvalues
        self.experimental_space_group = experimental_space_group
        self.experimental_unit_cell = experimental_unit_cell
        self.filtered_reflections = reflection_table
        self.filtered_reflections['z_value'] = zvalues
        
    def filter_data(self,variable_key,lower_limit_of_exclusion,upper_limit_of_exclusion):
        bad_data_selection=select_variable_range(self.filtered_reflections[variable_key],-1.0,0.0)
        inverse_sel = ~bad_data_selection
        self.filtered_reflections = self.filtered_reflections.select(inverse_sel)

    def map_indices_to_asu(self):
        crystal_symmetry=crystal.symmetry(unit_cell = self.experimental_unit_cell, space_group = self.experimental_space_group)
        miller_set=miller.set(crystal_symmetry = crystal_symmetry, indices = self.filtered_reflections['miller_index'])
        self.filtered_reflections["asu_miller_index"] = (miller_set.map_to_asu()).indices()



def select_variable_range(variable, lower_limit, upper_limit):
	sel = flex.bool()
	variable_aslist = list(variable)
	for j in range(0,len(variable_aslist)):
		if lower_limit < variable[j] <= upper_limit:
			sel.append(True)
		else:
			sel.append(False)
	Chosen_reflections = variable.select(sel)
	Chosen_indices = sel.iselection()
	#return Chosen_indices
	return sel