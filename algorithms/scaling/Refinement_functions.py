from dials_array_family_flex_ext import *
from cctbx.array_family.flex import *
from cctbx.array_family import flex
from scitbx import lbfgs

class optimiser:
    '''input data in form [Intensities,Variances,gvalues,sorted_l_index,sorted_h_index],[h_index_counter_array,h_index_cumulative_array]'''
    def __init__(self,data,h_index_arrays):
    	self.Intensities = data[0]
    	self.Variances = data[1]
    	self.x = data[2]
        self.sorted_l_index=data[3]
        self.sorted_h_index=data[4]
        self.h_index_arrays=h_index_arrays
        self.Ih_array = calc_Ih(self.x,self.Intensities,self.Variances,self.sorted_l_index,self.sorted_h_index,self.h_index_arrays)
        lbfgs.run(target_evaluator=self)

    def compute_functional_and_gradients(self):
        #self.Ih_array = calc_Ih(self.x,self.Intensities,self.Variances,self.sorted_l_index,self.sorted_h_index,self.h_index_arrays)
        self.Ih_array = calc_Ih(self.x,self.Intensities,self.Variances,self.sorted_l_index,self.sorted_h_index,self.h_index_arrays)
    	f = residual(self.Intensities,self.Variances,self.Ih_array,self.x,self.sorted_l_index,self.sorted_h_index)
    	g = gradient(self.Intensities,self.Variances,self.Ih_array,self.x,self.sorted_l_index,self.sorted_h_index,self.h_index_arrays)
        #print list(g)[0:10]
        #print list(self.x)[0:10]
        return f,g

	def callback_after_step(self,minimizer):
		print "LBFGS step"

def residual(Intensities,Variances,Ih_array,gvalues,sorted_l_index,sorted_h_index):
    R = 0.0
    for n in range(0,len(Intensities)):
        l=sorted_l_index[n]
        h=sorted_h_index[n]
        #if Variances[n]>0.0:
        R+=(((Intensities[n]-(gvalues[l-1]*Ih_array[h-1]))**2)/Variances[n])
        #else:
        #    R+=0.0
    print "R = " + str(R) 
    return R

def gradient(Intensities,Variances,Ih_array,gvalues,sorted_l_index,sorted_h_index,h_index_arrays):
    h_index_counter_array=h_index_arrays[0]
    h_index_cumulative_array=h_index_arrays[1]
    G = flex.double([0.0]*len(gvalues))
    dIhdg=calc_dIhdg(gvalues,Intensities,Variances,Ih_array,sorted_l_index,sorted_h_index,h_index_arrays)
    for n in range(0,len(h_index_counter_array)):
        lsum=h_index_counter_array[n]
        for i in range(0,lsum):
            #determine index of data table
            indexer=i+h_index_cumulative_array[n]
            l=sorted_l_index[indexer]
            h=sorted_h_index[indexer]
            #print l,h
            if l != 0 and h != 0:
                rhl = (Intensities[indexer]-(gvalues[l-1]*Ih_array[h-1]))/(Variances[indexer]**0.5)
                G[l-1] += 2.0*rhl*((-1.0*Ih_array[h-1]/(Variances[indexer]**0.5)) -((gvalues[l-1]/(Variances[indexer]**0.5))*dIhdg[l-1,h-1]))
            else:
                print l,h
    return G

def calc_dIhdg(gvalues,sorted_intensities,sorted_variances,Ih_array,sorted_l_index,sorted_h_index,h_index_arrays):
    h_index_counter_array,h_index_cumulative_array=h_index_arrays[0],h_index_arrays[1]
    dIh_array=flex.double([0.0]*len(gvalues)*len(h_index_counter_array))
    dIh_array.reshape(flex.grid(len(gvalues),len(h_index_counter_array)))
    numerator=flex.double([0.0]*len(gvalues)*len(h_index_counter_array))
    numerator.reshape(flex.grid(len(gvalues),len(h_index_counter_array)))
    denominator=flex.double([0.0]*len(h_index_counter_array))
    
    for n in range(0,len(h_index_counter_array)):
        a1=0.0
        b1=0.0
        lsum=h_index_counter_array[n]
        for i in range(0,lsum):
            #determine index of data table
            indexer=i+h_index_cumulative_array[n]
            l=sorted_l_index[indexer]
            h=sorted_h_index[indexer]
            if l != 0 and h != 0:
                #print l,h
            #print (sorted_variances[indexer], sorted_variances[indexer]**0.5)
                numerator[l-1,h-1] += ((sorted_intensities[indexer]/(sorted_variances[indexer])-(2.0*gvalues[l-1]*Ih_array[h-1]/(sorted_variances[indexer]))))
                denominator[h-1] += ((gvalues[l-1]**2)/sorted_variances[indexer])
    #print list(numerator)
    
    for h in range(0,len(h_index_counter_array)):
        if denominator[h] != 0.0:
            for l in range(0,len(gvalues)):
                dIh_array[l,h] = numerator[l,h]/denominator[h]
        else:
            print "zero denominator"
        #dIh_array.append(a1/b1)
        #else:
        #    dIh_array.append(0.0)
    return dIh_array

def calc_Ih(gvalues,sorted_intensities,sorted_variances,sorted_l_index,sorted_h_index,h_index_arrays):
	h_index_counter_array,h_index_cumulative_array=h_index_arrays[0],h_index_arrays[1]
	Ih_array=[]
	for h in range(0,len(h_index_counter_array)):
		a1=0.0
		b1=0.0
		lsum=h_index_counter_array[h]
		for i in range(0,lsum):
			#determine index of data table
			indexer=i+h_index_cumulative_array[h]
			#if sorted_variances[indexer] and sorted_intensities[indexer] > 0.0:
			l=sorted_l_index[indexer]
			a1+=(gvalues[l-1]*sorted_intensities[indexer]/sorted_variances[indexer])
			b1+=((gvalues[l-1]**2)/sorted_variances[indexer])
		#if b1!=0:
		Ih_array.append(a1/b1)
		#else:
		#	Ih_array.append(0.0)
	return Ih_array



