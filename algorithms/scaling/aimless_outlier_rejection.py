from dials.array_family import flex

def _reject_outliers(self, max_deviation):
  h_index_cumulative_array = self.Ih_table.h_index_cumulative_array
  outlier_list_h_index = []
  outlier_list_refl_index = []
  for i, n in enumerate(self.Ih_table.h_index_counter_array):
    if n > 2:
      outlier_found = first_test_for_an_outlier(self.Ih_table, i, n, max_deviation)
      if outlier_found:
        outlier_list_h_index.append(outlier_found[0])
        outlier_list_refl_index.append(outlier_found[1])
  if not outlier_list_h_index:
    return outlier_list_h_index, outlier_list_refl_index
  outlier_list_h_index, outlier_list_refl_index = iterative_test_for_subsequent_outliers(
    outlier_list_h_index, outlier_list_refl_index,
    h_index_cumulative_array, self.Ih_table, max_deviation)
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
  if n_new > 0:
    return iterative_test_for_subsequent_outliers(outlier_list_h_index,
      outlier_list_refl_index, h_index_cumulative_array, Ih_table, max_deviation)
  return outlier_list_h_index, outlier_list_refl_index

def subsequent_test_for_an_outlier(outlier_list_h_index, outlier_list_refl_index,
                                   h_index_cumulative_array, Ih_table, max_deviation):
  new_outliers_h_index = []
  new_outliers_refl_index = []
  Ihl = Ih_table.intensities
  scale_factors = Ih_table.inverse_scale_factors
  weights = Ih_table.weights
  for h_count_idx in outlier_list_h_index:
    h_idx_cumul = h_index_cumulative_array[h_count_idx:h_count_idx+2]
    refl_indices = range(h_idx_cumul[0], h_idx_cumul[1])
    for refl_idx in outlier_list_refl_index:
      if h_idx_cumul[0] <= refl_idx and refl_idx < h_idx_cumul[1]:
        refl_indices.remove(refl_idx)
    if len(refl_indices) > 2:
      norm_dev_list = flex.double([])
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



def first_test_for_an_outlier(Ih_table, i, n, max_deviation):
  h_idx_cumul = Ih_table.h_index_cumulative_array[i:i+2]
  I = Ih_table.intensities[h_idx_cumul[0]:h_idx_cumul[1]]
  g = Ih_table.inverse_scale_factors[h_idx_cumul[0]:h_idx_cumul[1]]
  w = Ih_table.weights[h_idx_cumul[0]:h_idx_cumul[1]]
  wgI = I*g*w
  wg2 = w*g*g
  fwd_prod = [0.0]
  rev_prod = [0.0]
  for j in range(n-1):
    fwd_prod.append(wgI[j] + fwd_prod[j])
    a = wgI[j] + rev_prod[0]
    rev_prod.insert(0, a)
  wgIsum = flex.double(fwd_prod) + flex.double(rev_prod)
  fwd_prod = [0.0]
  rev_prod = [0.0]
  for j in range(n-1):
    fwd_prod.append(wg2[j] + fwd_prod[j])
    a = wg2[j] + rev_prod[0]
    rev_prod.insert(0, a)
  wg2sum = flex.double(fwd_prod) + flex.double(rev_prod)
  norm_dev_list = (I - (g* wgIsum/wg2sum))/((1.0/w)+((g/wg2sum)**2))

  abs_norm_dev_list = (norm_dev_list**2)**0.5
  max_delta = abs_norm_dev_list.min_max_mean().max

  if max_delta > max_deviation:
    sel = norm_dev_list >= 0.0
    if sel.count(True) == 1:
      return (i, (sel.iselection() + Ih_table.h_index_cumulative_array[i])[0])
    elif sel.count(False) == 1:
      inv_sel = ~sel
      return (i, (inv_sel.iselection() + Ih_table.h_index_cumulative_array[i])[0])
    sel1 = abs(abs_norm_dev_list - max_delta) < 0.000001
    return (i, (sel1.iselection() + Ih_table.h_index_cumulative_array[i])[0])
  return None
