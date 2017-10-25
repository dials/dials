
from math import log10, floor

def calc_entropy(data, bin_boundaries):
  entropy=0.0
  for i in range(len(bin_boundaries)-1):
    c=0.0
    lowerbound = bin_boundaries[i]
    upperbound = bin_boundaries[i+1]
    for j in data:
      if lowerbound <= j < upperbound:
        c=c+1.0
    if c>0.0:
      entropy += (c/float(len(data)))*log10(c/float(len(data)))
  print entropy

def calc_entropy_2D(refl_table, x_bin_boundaries, y_bin_boundaries):
  entropy=0.0
  #loop through detector areas
  for i in range(len(x_bin_boundaries)-1):
    x_lowerbound = x_bin_boundaries[i]
    x_upperbound = x_bin_boundaries[i+1]
    for j in range(len(y_bin_boundaries)-1):
      y_lowerbound = y_bin_boundaries[j]
      y_upperbound = y_bin_boundaries[j+1]
      c_array_in_z = [0]*refl_table.nzbins
      for n in range(len(refl_table.h_index_counter)):
        #loop over groups of equivalent reflections
        c=0.0
        data_x = refl_table.dataset_x[refl_table.h_index_cumulative[n]:refl_table.h_index_cumulative[n+1]]
        data_y = refl_table.dataset_y[refl_table.h_index_cumulative[n]:refl_table.h_index_cumulative[n+1]]
        for k in range(len(data_x)):
          #calculate entropy due to member of each group
          if x_lowerbound <= data_x[k] < x_upperbound:
            if y_lowerbound <= data_y[k] < y_upperbound:
              index = refl_table.h_index_cumulative[n] + k
              z_bin = refl_table.z_bins[index]
              c_array_in_z[z_bin] += 1.0
              #c=c+1.0
        for c in c_array_in_z:
          if c>0.0:
            entropy += (c/float(len(data_x)))*log10(c/float(len(data_x)))
  return entropy

class refl_table(object):
  def __init__(self, dataset_x, dataset_y, z_bins, nzbins, h_index_counter, h_index_cumulative):
    self.dataset_x = dataset_x
    self.dataset_y = dataset_y
    self.z_bins = z_bins
    self.nzbins = nzbins
    self.h_index_counter = h_index_counter
    self.h_index_cumulative = h_index_cumulative

def count_members_in_child(reflection_table, bin_object):
  count = [0.0]*reflection_table.nzbins
  for k in range(len(reflection_table.dataset_x)):
    if bin_object.x_bin_bounds[0] <= reflection_table.dataset_x[k] < bin_object.x_bin_bounds[1]:
      if bin_object.y_bin_bounds[0] <= reflection_table.dataset_y[k] < bin_object.y_bin_bounds[1]:
        z_bin = reflection_table.z_bins[k]
        count[z_bin] += 1.0
  return count

def try_a_divide(reflection_table, bin_obj, axis, min_in_each_bin):
  starting_entropy = calc_entropy_2D(reflection_table, bin_obj.x_bin_bounds, bin_obj.y_bin_bounds)
  if axis == 'x':
    bin_obj.new_x_bounds = bin_obj.x_bin_bounds[:]
    bin_obj.new_x_bounds.insert(1, bin_obj.x_bin_bounds[0] + 1)
    new_entropy = calc_entropy_2D(reflection_table, bin_obj.new_x_bounds, bin_obj.y_bin_bounds)
    if new_entropy < starting_entropy:
      bin_obj.children.append(detector_bin_child(bin_obj.number, bin_obj.new_x_bounds[1:], bin_obj.y_bin_bounds))
      bin_obj.children.append(detector_bin_child(bin_obj.number, bin_obj.new_x_bounds[0:2], bin_obj.y_bin_bounds))
      members_in_child_1 = count_members_in_child(reflection_table, bin_obj.children[0])
      members_in_child_2 = count_members_in_child(reflection_table, bin_obj.children[1])
      if min(members_in_child_1) < min_in_each_bin or min(members_in_child_2) < min_in_each_bin:
          bin_obj.children = []
      else:
        print "good divide in %s for bin %s" % (axis, bin_obj.number)
  elif axis == 'y':
    bin_obj.new_y_bounds = bin_obj.y_bin_bounds[:]
    bin_obj.new_y_bounds.insert(1, bin_obj.y_bin_bounds[0] + 1)
    new_entropy = calc_entropy_2D(reflection_table, bin_obj.x_bin_bounds, bin_obj.new_y_bounds)
    if new_entropy < starting_entropy:
      bin_obj.children.append(detector_bin_child(bin_obj.number, bin_obj.x_bin_bounds, bin_obj.new_y_bounds[1:]))
      bin_obj.children.append(detector_bin_child(bin_obj.number, bin_obj.x_bin_bounds, bin_obj.new_y_bounds[0:2]))
      members_in_child_1 = count_members_in_child(reflection_table, bin_obj.children[0])
      members_in_child_2 = count_members_in_child(reflection_table, bin_obj.children[1])
      if min(members_in_child_1) < min_in_each_bin or min(members_in_child_2) < min_in_each_bin:
          bin_obj.children = []
      else:
        print "good divide in %s for bin %s" % (axis, bin_obj.number)



class detector_bin(object):
  def __init__(self, number, bin_boundaries_x1, bin_boundaries_y1):
    self.number = number
    self.x_coord = self.number%3
    self.y_coord = int(floor(self.number/3))
    self.x_bin_bounds = [bin_boundaries_x1[self.x_coord],bin_boundaries_x1[self.x_coord+1]]
    self.y_bin_bounds = [bin_boundaries_y1[self.y_coord],bin_boundaries_y1[self.y_coord+1]]
    self.children = []

class detector_bin_child(object):
  def __init__(self, number, bin_boundaries_x1, bin_boundaries_y1):
    self.number = number
    self.x_bin_bounds = bin_boundaries_x1
    self.y_bin_bounds = bin_boundaries_y1
    self.children = []


def perform_optimal_divide(reflection_table, bin_boundaries_x1, bin_boundaries_y1, min_in_each_bin):
  binlist=[]
  for i in range(0,9):
    bin_object = detector_bin(i, bin_boundaries_x1, bin_boundaries_y1)
    try_a_divide(reflection_table, bin_object, 'y', min_in_each_bin)
    binlist.append(bin_object)

  for bin in binlist:
    if len(bin.children)==0:
      try_a_divide(reflection_table, bin, 'x', min_in_each_bin)
    else:
      for child in bin.children:
        try_a_divide(reflection_table, child, 'x', min_in_each_bin)

  final_bin_grid = []

  for bin in binlist:
    if len(bin.children)==0:
      print "large bin %s" % bin.number
      final_bin_grid.append((bin.x_bin_bounds,bin.y_bin_bounds))
    else:
      for child in bin.children:
        if len(child.children)==0:
          print "child bin %s" % child.number
          final_bin_grid.append((child.x_bin_bounds,child.y_bin_bounds))
        else:
          for grandchild in child.children:
            print "grandchild bin %s" % grandchild.number
            final_bin_grid.append((grandchild.x_bin_bounds,grandchild.y_bin_bounds))

  return final_bin_grid


if __name__ == "__main__":

  #simulate two reflection groups
  x_data = [2, 0, 1, 1, 2, 2, 3, 3, 3, 4, 4]
  y_data = [1, 3, 3, 2, 3, 3, 2, 2, 4, 2, 3]

  z_bins = [0,0,0,0,0,0,0,0,0,0,0,0,1,1]

  x_data_2 = [1, 4, 5]
  y_data_2 = [1, 3, 3]

  dataset_x = x_data + x_data_2
  dataset_y = y_data + y_data_2
  h_index_counter = [len(x_data), len(x_data_2)]
  h_index_cumulative = [0, h_index_counter[0], h_index_counter[0] + h_index_counter[1]]

  #build a data_manger object with reflection table-like organisation
  reflection_table = refl_table(dataset_x, dataset_y, z_bins, 2, h_index_counter, h_index_cumulative)

  bin_boundaries_x1 = [0,2,4,6]
  bin_boundaries_y1 = [0,2,4,6]

  print calc_entropy_2D(reflection_table, bin_boundaries_x1, bin_boundaries_y1)

  optimal_boundaries = perform_optimal_divide(reflection_table, bin_boundaries_x1, bin_boundaries_y1, min_in_each_bin=1)
  print optimal_boundaries

  reflection_table.a_bin_index = [-1 for i in reflection_table.dataset_x]


  for i, x_reflection in enumerate(reflection_table.dataset_x):
    y_reflection = reflection_table.dataset_y[i]
    for j, boundaries in enumerate(optimal_boundaries):
      if boundaries[0][0] <= x_reflection < boundaries[0][1] and boundaries[1][0] <= y_reflection < boundaries[1][1]:
        reflection_table.a_bin_index[i] = j
  print reflection_table.a_bin_index
