import pickle as pickle
from dials_array_family_flex_ext import *
from cctbx.array_family.flex import *
from cctbx.array_family import flex
from cctbx import miller
from cctbx import crystal
import matplotlib.pyplot as plt
import numpy as np
from dials.util.options import Importer, flatten_experiments, flatten_reflections, OptionParser


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


'''load experiments.json file and reflection table'''
filepath="/Users/whi10850/Documents/test_data/integrate/13_integrated.pickle"
json_filepath="/Users/whi10850/Documents/test_data/integrate/13_integrated_experiments.json"
parsestring = OptionParser(read_experiments=True,read_reflections=True,check_format=False)
params, options = parsestring.parse_args([json_filepath,filepath])
experiments = flatten_experiments(params.input.experiments)
reflections_table = flatten_reflections(params.input.reflections)[0]
experimental_unit_cell=experiments.crystals()[0].get_unit_cell().parameters()
experimental_space_group=experiments.crystals()[0].get_space_group()


'''extract dvalues, miller indices and zvalues'''
dvalues=reflections_table['d']
miller_indices=reflections_table['miller_index']
xyzvalues=reflections_table['xyzobs.px.value']
zvalues=flex.double([xyzvalues[i][2] for i in range(0,len(xyzvalues))])
Ivalues=reflections_table['intensity.prf.value']
sigmavalues=reflections_table['intensity.prf.variance']


'''determine miller indices in asymmetric unit'''
crystal_symmetry=crystal.symmetry(
    unit_cell=experimental_unit_cell,
	space_group=experimental_space_group)
miller_set=miller.set(crystal_symmetry=crystal_symmetry,
	indices=reflections_table['miller_index'])
asu_miller_indices=miller_set.map_to_asu()

'''determine range for gridding of the 'decay' correction'''
zmax=max(flex.double(zvalues))
zmin=min(flex.double(zvalues))
dmax=max(flex.double(dvalues))
dmin=min(flex.double(dvalues))
print 'dmin,dmax,zmin,zmax = '+'%f, %f, %f,   %f' % (dmin,dmax,zmin,zmax)

ndbins=20
nzbins=20
nlbins=ndbins*nzbins

'''attempt to creating equally populated d bins - values determined here'''
def calculate_evenly_populated_resolution_bins(dvalues,asu_miller_indices):
	permuted = asu_miller_indices.sort_permutation(by_value='resolution',reverse=True)
	number_of_entries=len(dvalues)
	dvalues_permuted=dvalues.select(permuted)
	resolution_bins=flex.double([0]*(ndbins+1))
	for i in range(0,ndbins):
		n=(i*number_of_entries//ndbins)
		resolution_bins[i]=dvalues_permuted[n]
	resolution_bins[ndbins]=dmax
	return resolution_bins

resolution_bins=calculate_evenly_populated_resolution_bins(dvalues,asu_miller_indices)
#resolution_bins=flex.double(range(0,21))*dmax/20.0
z_bins=flex.double(range(0,nzbins+1))*zmax/nzbins

'''calculate index for the gridding'''
d_bin_index=flex.int([0]*len(dvalues))
z_bin_index=flex.int([0]*len(zvalues))
for i in range(0,ndbins):
	dselection=select_variable_range(dvalues,resolution_bins[i],resolution_bins[i+1])
	d_bin_index.set_selected(dselection,i+1)
for i in range(0,nzbins):
	zselection=select_variable_range(zvalues,z_bins[i],z_bins[i+1])
	z_bin_index.set_selected(zselection,i+1)
dz_bin_index=(d_bin_index+((z_bin_index-1)*ndbins))

'''reorder columns of interest so that they can be iterated over'''
permuted = asu_miller_indices.sort_permutation(by_value='packed_indices')
sorted_miller_indices=(asu_miller_indices.indices()).select(permuted)
sorted_h_index=flex.int([0]*len(permuted))
sorted_h_index[0]=1

'''build the h_index array and also the h_index counter array i.e. how many reflections for each unique reflection'''
h_index_counter_array=[]
h_index=1
h_index_counter=1
for i in range(1,len(permuted)):
	if sorted_miller_indices[i]==sorted_miller_indices[i-1]:
		sorted_h_index[i]=h_index
		h_index_counter+=1
	else:
		h_index+=1
		sorted_h_index[i]=h_index
		h_index_counter_array.append(h_index_counter)
		h_index_counter=1
h_index_counter_array.append(h_index_counter)

'''calculate the cumulative sum after each h_index group'''
hsum=0
h_index_cumulative_array=[0]
for h in range(0,len(h_index_counter_array)):
	hsum+=h_index_counter_array[h]
	h_index_cumulative_array.append(hsum)

sorted_l_index=dz_bin_index.select(permuted)
sorted_intensities=Ivalues.select(permuted)
sorted_variances=sigmavalues.select(permuted)

'''build data and index arrays to pass to functions for minimisation calculations'''
sorted_data=[sorted_miller_indices,sorted_intensities,sorted_variances,sorted_l_index,sorted_h_index]
h_index_arrays=[h_index_counter_array,h_index_cumulative_array]

G0_values=flex.double([1.0]*(ndbins*nzbins))


def calc_Ih(gl,G_values,sorted_data,h_index_arrays):
	sorted_miller_indices,sorted_intensities,sorted_variances,sorted_l_index,sorted_h_index=sorted_data[0],sorted_data[1],sorted_data[2],sorted_data[3],sorted_data[4]
	h_index_counter_array,h_index_cumulative_array=h_index_arrays[0],h_index_arrays[1]
	Ih_array=[]
	for h in range(0,len(h_index_counter_array)):
		a1=0.0
		b1=0.0
		lsum=h_index_counter_array[h]
		for i in range(0,lsum):
			#determine index of data table
			indexer=i+h_index_cumulative_array[h]
			if sorted_variances[indexer] and sorted_intensities[indexer] > 0.0:
				l=sorted_l_index[indexer]
				a1+=(gl[l-1]*sorted_intensities[indexer]*G_values[l-1]/sorted_variances[indexer])
				b1+=(((gl[l-1]*G_values[l-1])**2)/sorted_variances[indexer])
		if b1!=0:
			Ih_array.append(a1/b1)
		else:
			Ih_array.append(0.0)
	return Ih_array

sigma=0.05

def update_gl(G_values, Ih_array, sorted_data, h_index_arrays):
	sorted_miller_indices,sorted_intensities,sorted_variances,sorted_l_index,sorted_h_index=sorted_data[0],sorted_data[1],sorted_data[2],sorted_data[3],sorted_data[4]
	h_index_counter_array,h_index_cumulative_array=h_index_arrays[0],h_index_arrays[1]
	gl_array=flex.double([0.0]*nlbins)
	gl_numerator=flex.double([0.0]*nlbins)
	gl_denominator=flex.double([0.0]*nlbins)
	R = residual(sorted_intensities,sorted_variances,Ih_array,G_values,sorted_l_index,sorted_h_index)
	for n in range(0,len(sorted_intensities)):
		l=sorted_l_index[n]
		if sorted_variances[n] and sorted_intensities[n] > 0.0:
			gl_numerator[l-1]+=(Ih_array[sorted_h_index[n]-1]*sorted_intensities[n]*G_values[l-1]/sorted_variances[n])
			gl_denominator[l-1]+=((Ih_array[sorted_h_index[n]-1]**2)*(G_values[l-1]**2)/sorted_variances[n])
	for l in range(0,len(gl_array)):
		gl_array[l-1]=((gl_numerator[l-1]+(G_values[l-1]/(sigma**2)))/(gl_denominator[l-1]+((G_values[l-1]/sigma)**2)))
	return gl_array

gl0_array=flex.double([1.0]*nlbins)
Ih0_array=calc_Ih(gl0_array,G0_values,sorted_data,h_index_arrays)

def residual(Intensities,Variances,Ih_array,gvalues,sorted_l_index,sorted_h_index):
    R = 0.0
    for n in range(0,len(Intensities)):
        l=sorted_l_index[n]
        h=sorted_h_index[n]
        if Variances[n]>0.0:
        	R+=(((Intensities[n]-(gvalues[l-1]*Ih_array[h-1]))**2)/Variances[n])
        #else:
        #    R+=0.0
    print "R = " + str(R) 
    return R

def update_G_values(Gl_array,Ih_array,counter,niter,sorted_data,h_index_arrays):
	gl_new = update_gl(Gl_array,Ih_array,sorted_data,h_index_arrays)
	Ih_new = calc_Ih(gl_new,Gl_array,sorted_data,h_index_arrays)
	G_new = gl_new*Gl_array
	counter+=1
	if counter<niter:
		return update_G_values(G_new,Ih_new,counter,niter,sorted_data,h_index_arrays)
	else:
		return G_new, Ih_new
		

'''do a minimisation procedure'''		
G_fin, Ih_fin = update_G_values(G0_values,Ih0_array,0,20,sorted_data,h_index_arrays)


'''print out various values for inspection'''
print "Miller Index, l_index, G_l, Measured Intensity, Measured Variance, ""True"" Intensity"
for i in range(15140,15150):
    print (list(sorted_miller_indices)[i],sorted_l_index[i],list(G_fin)[sorted_l_index[i]-1],list(sorted_intensities)[i],list(sorted_variances)[i],list(Ih_fin)[sorted_h_index[i]-1])


G_fin_list=list(G_fin.as_double())
print len(G_fin_list)
G_fin_2d=np.reshape(G_fin_list,(nzbins,ndbins)).T
im=plt.imshow(G_fin_2d,cmap='viridis',origin='lower')
plt.colorbar(im)
plt.ylabel('resolution (d) bin')
plt.xlabel('time (z) bin')
plt.title('$G_l$ correction factors using Kabsch method')
plt.savefig('Scaling_output_figure.png')
plt.show()


plt.figure(1)
plt.plot(list(G_fin)[0:20])
plt.plot(list(G_fin)[20:40])
plt.figure(2)
plt.plot(list(G_fin)[2:382:20])
plt.plot(list(G_fin)[8:388:20])
plt.show()



exit()
