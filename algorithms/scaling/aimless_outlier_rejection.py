from dials.array_family import flex

def reject_outliers(self, max_deviation):
  Ih_table = self.Ih_table.Ih_table
  h_index_cumulative_array = self.Ih_table.h_index_cumulative_array
  outlier_list_h_index=[]
  outlier_list_refl_index=[]
  for i, n in enumerate(self.Ih_table.h_index_counter_array):
    #index = h_index_cumulative_array[i]
    if n > 2:
      h_idx_cumul = h_index_cumulative_array[i:i+2]
      Ihls_u = Ih_table['intensity'][h_idx_cumul[0]:h_idx_cumul[1]]
      gs_u = Ih_table['inverse_scale_factor'][h_idx_cumul[0]:h_idx_cumul[1]]
      ws_u = Ih_table['weights'][h_idx_cumul[0]:h_idx_cumul[1]]
      outlier_found = first_test_for_an_outlier(h_index_cumulative_array, Ihls_u, 
        gs_u, ws_u, i, n, max_deviation)
      if outlier_found:
        outlier_list_h_index.append(outlier_found[0])
        outlier_list_refl_index.append(outlier_found[1])
  if not outlier_list_h_index:
    return outlier_list_h_index, outlier_list_refl_index
  else:
    outlier_list_h_index, outlier_list_refl_index = iterative_test_for_subsequent_outliers(
      outlier_list_h_index, outlier_list_refl_index,
      h_index_cumulative_array, Ih_table, max_deviation)
    return outlier_list_h_index, outlier_list_refl_index


def iterative_test_for_subsequent_outliers(outlier_list_h_index, outlier_list_refl_index,
    h_index_cumulative_array, Ih_table, max_deviation):
  n_found = len(outlier_list_h_index)
  new_outliers_h_index, new_outliers_refl_index = subsequent_test_for_an_outlier(
    outlier_list_h_index, outlier_list_refl_index,
    h_index_cumulative_array, Ih_table, max_deviation)
  for new_outlier, new_outliers_h in zip(new_outliers_refl_index, new_outliers_h_index):
    if new_outlier not in outlier_list_refl_index:
      outlier_list_h_index.append(new_outliers_h)
      outlier_list_refl_index.append(new_outlier)
  n_new = len(outlier_list_h_index) - n_found
  if n_new>0:
    #outlier_list_h_index = outlier_list_h_index + new_outliers_h_index
    #outlier_list_refl_index = outlier_list_refl_index + new_outliers_refl_index
    return iterative_test_for_subsequent_outliers(outlier_list_h_index, outlier_list_refl_index,
                                   h_index_cumulative_array, Ih_table, max_deviation)
  else:
    return outlier_list_h_index, outlier_list_refl_index

def subsequent_test_for_an_outlier(outlier_list_h_index, outlier_list_refl_index,
                                   h_index_cumulative_array, Ih_table, max_deviation):
  new_outliers_h_index=[]
  new_outliers_refl_index=[]
  Ihl = Ih_table['intensity']
  scale_factors = Ih_table['inverse_scale_factor']
  weights = Ih_table['weights']
  #print outlier_list_h_index
  #print outlier_list_refl_index                                 
  for h_count_idx in outlier_list_h_index:
    h_idx_cumul = h_index_cumulative_array[h_count_idx:h_count_idx+2]
    refl_indices = range(h_idx_cumul[0],h_idx_cumul[1])
    #print refl_indices
    for refl_idx in outlier_list_refl_index:
      if h_idx_cumul[0] <= refl_idx and refl_idx < h_idx_cumul[1]:
        #print refl_idx
        refl_indices.remove(refl_idx)
    if len(refl_indices) > 2:
      norm_dev_list=flex.double([])
      for index in refl_indices:
        Is = flex.double([])
        ws = flex.double([])
        gs = flex.double([])
        for index_2 in refl_indices:
          if index_2 != index:
            Is.append(Ihl[index_2])
            ws.append(weights[index_2])
            gs.append(scale_factors[index_2])
        Ih = flex.sum(ws*gs*Is)/flex.sum(ws*(gs**2))
        sigma = 1.0/flex.sum(ws*(gs**2))
        norm_dev_list.append((Ihl[index] - scale_factors[index]*Ih)/
          (1.0/weights[index] + ((scale_factors[index]*sigma)**2)))
      abs_norm_dev_list = abs(norm_dev_list)
      max_delta = max(abs_norm_dev_list)
      if max_delta > max_deviation:
        #there is an outlier, so proceed to find it and rerun search
        sel = norm_dev_list >= 0.0
        n_pos = sel.count(True)
        n_neg = sel.count(False)
        if n_pos == 1:
          m = sel.iselection()
          idx = m + h_index_cumulative_array[h_count_idx]
          new_outliers_h_index.append(h_count_idx)
          new_outliers_refl_index.append(idx[0])
        elif n_neg == 1:
          inv_sel = ~sel
          m = inv_sel.iselection()
          idx = m + h_index_cumulative_array[h_count_idx]
          new_outliers_h_index.append(h_count_idx)
          new_outliers_refl_index.append(idx[0])
        else: 
          sel1 = abs(abs_norm_dev_list - max_delta) < 0.000001
          m = sel1.iselection()
          idx = m + h_index_cumulative_array[h_count_idx]
          new_outliers_h_index.append(h_count_idx)
          new_outliers_refl_index.append(idx[0])
  return new_outliers_h_index, new_outliers_refl_index



def first_test_for_an_outlier(h_index_cumulative_array, Ihls_u, gs_u, ws_u, i, n, max_deviation):
  norm_dev_list = flex.double([])
  #print ws_u.count(0.0)
  for j in range(n): 
    Is = Ihls_u[0:j]
    ws = ws_u[0:j]
    gs = gs_u[0:j]
    Is.extend(Ihls_u[j+1:n])
    ws.extend(ws_u[j+1:n])
    gs.extend(gs_u[j+1:n])
    sigma = 1.0/flex.sum(ws*(gs**2))
    Ih = flex.sum(ws*gs*Is) * sigma
    norm_dev_list.append((Ihls_u[j] - (gs_u[j]*Ih))/
                         ((1.0/ws_u[j]) + ((gs_u[j]*sigma)**2)))
  abs_norm_dev_list = abs(norm_dev_list)
  max_delta = max(abs_norm_dev_list)
  if max_delta > max_deviation:
    sel = norm_dev_list >= 0.0
    n_pos = sel.count(True)
    n_neg = sel.count(False)
    if n_pos == 1:
      m = sel.iselection()
      idx = m + h_index_cumulative_array[i]
      return (i,idx[0])
    elif n_neg == 1:
      inv_sel = ~sel
      m = inv_sel.iselection()
      idx = m + h_index_cumulative_array[i]
      return (i,idx[0])
    else: 
      sel1 = abs(abs_norm_dev_list - max_delta) < 0.000001
      m = sel1.iselection()
      idx = m + h_index_cumulative_array[i]
      return (i,idx[0])
  else:
    return None