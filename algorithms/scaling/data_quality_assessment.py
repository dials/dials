def R_meas(Data_Manager_object):
    Ihl = Data_Manager_object.sorted_reflections['intensity.sum.value']
    h_array = Data_Manager_object.sorted_reflections['h_index']
    l_array = Data_Manager_object.sorted_reflections['l_bin_index']
    a_array = Data_Manager_object.sorted_reflections['a_bin_index']
    xy_array = Data_Manager_object.sorted_reflections['xy_bin_index']
    Ih_array = Data_Manager_object.Ih_array
    g_values = Data_Manager_object.g_values
    g2_values = Data_Manager_object.g2_values
    g3_values = Data_Manager_object.g3_values

    Rmeas_upper = 0.0
    Rmeas_lower = 0.0
    ndbins = Data_Manager_object.ndbins

    for h in range(Data_Manager_object.n_unique_indices):
        a1 = 0.0
        lsum = Data_Manager_object.h_index_counter_array[h]
        if lsum > 1:
            for i in range(lsum):
                indexer = i + Data_Manager_object.h_index_cumulative_array[h]
                (h,l,a,xy)=Data_Manager_object.bin_index[indexer]
                #l = l_array[indexer]
                #h = h_array[indexer]
                #a = a_array[indexer]
                #xy = xy_array[indexer]
                a1 += abs((Ihl[indexer]/(g_values[l]*g2_values[a]*g3_values[xy])) - Ih_array[h])
                Rmeas_lower += (Ihl[indexer]/(g_values[l]*g2_values[a]*g3_values[xy]))
            Rmeas_upper += (((float(lsum)/(float(lsum)-1.0))**0.5) * a1)
    
    Rmeas = Rmeas_upper/Rmeas_lower
    return Rmeas    

def R_pim(Data_Manager_object):
    Ihl = Data_Manager_object.sorted_reflections['intensity.sum.value']
    h_array = Data_Manager_object.sorted_reflections['h_index']
    l_array = Data_Manager_object.sorted_reflections['l_bin_index']
    a_array = Data_Manager_object.sorted_reflections['a_bin_index']
    xy_array = Data_Manager_object.sorted_reflections['xy_bin_index']
    Ih_array = Data_Manager_object.Ih_array
    g_values = Data_Manager_object.g_values
    g2_values = Data_Manager_object.g2_values
    g3_values = Data_Manager_object.g3_values
    Rpim_upper = 0.0
    Rpim_lower = 0.0
    ndbins = Data_Manager_object.ndbins
    for h in range(Data_Manager_object.n_unique_indices):
        a1 = 0.0
        lsum = Data_Manager_object.h_index_counter_array[h]
        if lsum > 1:
            for i in range(lsum):
                indexer = i + Data_Manager_object.h_index_cumulative_array[h]
                (h,l,a,xy)=Data_Manager_object.bin_index[indexer]
                #l = l_array[indexer]
                #h = h_array[indexer]
                #a = a_array[indexer]
                #xy = xy_array[indexer]
                a1 += abs((Ihl[indexer]/(g_values[l]*g2_values[a]*g3_values[xy])) - Ih_array[h])
                Rpim_lower += (Ihl[indexer]/(g_values[l]*g2_values[a]*g3_values[xy]))
            Rpim_upper += (((1.0/(float(lsum)-1.0))**0.5) * a1)

    Rpim = Rpim_upper/Rpim_lower
    return Rpim          


