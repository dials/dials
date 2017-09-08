'''
Simple functions for calculating R_meas, R_int from a Data_Manager_object
'''

def R_meas(Data_Manager_object):
    '''Calculate R_meas from a Data_Manager_object'''
    Ihl = Data_Manager_object.sorted_reflections[Data_Manager_object.int_method[0]]
    Ih = Data_Manager_object.Ih_values
    '''gvalues1 = flex.double([Data_Manager_object.g_values[i]
                for i in Data_Manager_object.sorted_reflections['l_bin_index']])
    gvalues2 = flex.double([Data_Manager_object.g2_values[i]
                for i in Data_Manager_object.sorted_reflections['a_bin_index']])
    gvalues3 = flex.double([Data_Manager_object.g3_values[i]
                for i in Data_Manager_object.sorted_reflections['xy_bin_index']])

    res = abs((Ihl/(gvalues1*gvalues2*gvalues3))- Ih)
    sumres = np.bincount(self.data_manager.sorted_reflections['h_index'], res)
    nh = Data_Manager_object.h_index_counter_array
    prefactor = (nh/(nh-1.0))**0.5
    gsq = (((self.gproduct)**2)/variances)
        sumgsq = np.bincount(self.data_manager.sorted_reflections['h_index'], gsq)'''

    g_values = Data_Manager_object.g_values
    g2_values = Data_Manager_object.g2_values
    g3_values = Data_Manager_object.g3_values

    Rmeas_upper = 0.0
    Rmeas_lower = 0.0
    for h in range(len(Data_Manager_object.h_index_counter_array)):
        a1 = 0.0
        lsum = Data_Manager_object.h_index_counter_array[h]
        if lsum > 1:
            for i in range(lsum):
                indexer = i + Data_Manager_object.h_index_cumulative_array[h]
                l = Data_Manager_object.sorted_reflections['l_bin_index'][indexer]
                a = Data_Manager_object.sorted_reflections['a_bin_index'][indexer]
                xy = Data_Manager_object.sorted_reflections['xy_bin_index'][indexer]
                a1 += abs((Ihl[indexer] / (g_values[l] * g2_values[a] * g3_values[xy])) - Ih[indexer])
                Rmeas_lower += (Ihl[indexer] / (g_values[l] * g2_values[a] * g3_values[xy]))
            Rmeas_upper += (((float(lsum) / (float(lsum) - 1.0))**0.5) * a1)
    Rmeas = Rmeas_upper / Rmeas_lower
    return Rmeas

def R_pim(Data_Manager_object):
    '''Calculate R_pim from a Data_Manager_object'''
    Ihl = Data_Manager_object.sorted_reflections[Data_Manager_object.int_method[0]]
    Ih = Data_Manager_object.Ih_values
    g_values = Data_Manager_object.g_values
    g2_values = Data_Manager_object.g2_values
    g3_values = Data_Manager_object.g3_values
    Rpim_upper = 0.0
    Rpim_lower = 0.0
    for h in range(len(Data_Manager_object.h_index_counter_array)):
        a1 = 0.0
        lsum = Data_Manager_object.h_index_counter_array[h]
        if lsum > 1:
            for i in range(lsum):
                indexer = i + Data_Manager_object.h_index_cumulative_array[h]
                l = Data_Manager_object.sorted_reflections['l_bin_index'][indexer]
                a = Data_Manager_object.sorted_reflections['a_bin_index'][indexer]
                xy = Data_Manager_object.sorted_reflections['xy_bin_index'][indexer]

                a1 += abs((Ihl[indexer] / (g_values[l] * g2_values[a] * g3_values[xy])) - Ih[indexer])
                Rpim_lower += (Ihl[indexer] / (g_values[l] * g2_values[a] * g3_values[xy]))
            Rpim_upper += (((1.0 / (float(lsum) - 1.0))**0.5) * a1)
    Rpim = Rpim_upper / Rpim_lower
    return Rpim
