def R_meas(Data_Manager_object):
    Ihl = Data_Manager_object.sorted_reflections['intensity.sum.value']
    h_array = Data_Manager_object.sorted_reflections['h_index']
    l_array = Data_Manager_object.sorted_reflections['l_bin_index']
    Ih_array = Data_Manager_object.Ih_array
    g_values = Data_Manager_object.g_values

    Rmeas_upper = 0.0
    Rmeas_lower = 0.0
    ndbins = Data_Manager_object.ndbins

    for h in range(Data_Manager_object.n_unique_indices):
        a1 = 0.0
        lsum = Data_Manager_object.h_index_counter_array[h]
        if lsum > 1:
            for i in range(lsum):
                indexer = i + Data_Manager_object.h_index_cumulative_array[h]
                l = l_array[indexer]
                h = h_array[indexer]
                a1 += abs((Ihl[indexer]/g_values[l]) - Ih_array[h])
                Rmeas_lower += Ihl[indexer]
            Rmeas_upper += (((float(lsum)/(float(lsum)-1.0))**0.5) * a1)
    
    Rmeas = Rmeas_upper/Rmeas_lower

    return Rmeas    

def R_pim(Data_Manager_object):
    Ihl = Data_Manager_object.sorted_reflections['intensity.sum.value']
    h_array = Data_Manager_object.sorted_reflections['h_index']
    l_array = Data_Manager_object.sorted_reflections['l_bin_index']
    Ih_array = Data_Manager_object.Ih_array
    g_values = Data_Manager_object.g_values
    Rpim_upper = 0.0
    Rpim_lower = 0.0
    ndbins = Data_Manager_object.ndbins
    for h in range(Data_Manager_object.n_unique_indices):
        a1 = 0.0
        lsum = Data_Manager_object.h_index_counter_array[h]
        if lsum > 1:
            for i in range(lsum):
                indexer = i + Data_Manager_object.h_index_cumulative_array[h]
                l = l_array[indexer]
                h = h_array[indexer]
                a1 += abs((Ihl[indexer]/g_values[l]) - Ih_array[h])
                Rpim_lower += Ihl[indexer]
            Rpim_upper += (((1.0/(float(lsum)-1.0))**0.5) * a1)
    
    Rpim = Rpim_upper/Rpim_lower
    return Rpim          


