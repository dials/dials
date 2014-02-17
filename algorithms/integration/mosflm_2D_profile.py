from __future__ import division
from dials.model.data import Reflection, ReflectionList
from dials.algorithms.integration import add_2d, subtrac_bkg_2d, \
              fitting_2d_partials, fitting_2d_multile_var_build_mat, sigma_2d


from scitbx.array_family import flex

def make_2d_profile(reflections):
  #print "len(reflections) =", len(reflections)

  max_i_01 = 0.0
  for ref in reflections:
    if ref.is_valid():
      if ref.intensity > max_i_01:
        max_i_01 = ref.intensity
  #print "max_i_01 =", max_i_01
  max_i = 0.0
  for ref in reflections:
    if ref.is_valid():
      if ref.intensity > max_i and ref.intensity < max_i_01 * 0.95:
        max_i = ref.intensity
  thold = 0.5 * max_i

  select_rlist = ReflectionList()
  for ref in reflections:
    if ref.is_valid() and ref.intensity > thold and ref.intensity < max_i:
      select_rlist.append(ref)
  counter = 0
  #print "len(select_rlist) =", len(select_rlist)
  big_nrow = 0
  big_ncol = 0
  for ref in select_rlist:
    local_nrow = ref.shoebox.all()[1]
    local_ncol = ref.shoebox.all()[2]
    if local_nrow > big_nrow:
      big_nrow = local_nrow
    if local_ncol > big_ncol:
      big_ncol = local_ncol
    counter += 1
  #print big_nrow, big_ncol
  big_nrow = big_nrow * 2 + 1
  big_ncol = big_ncol * 2 + 1
  sumation = flex.double(flex.grid(big_nrow, big_ncol), 0)
  descr = flex.double(flex.grid(1, 3), 0)
  for ref in select_rlist:
    shoebox = ref.shoebox
    #mask = ref.shoebox_mask                                 # may be needed soon
    background = ref.shoebox_background
    data2d = shoebox[0:1, :, :]
    #mask2d = mask[0:1, :, :]                                # may be needed soon
    background2d = background[0:1, :, :]
    data2d.reshape(flex.grid(shoebox.all()[1:]))
    #mask2d.reshape(flex.grid(shoebox.all()[1:]))            # may be needed soon
    background2d.reshape(flex.grid(shoebox.all()[1:]))

    descr[0, 0] = ref.centroid_position[0] - ref.bounding_box[0]
    descr[0, 1] = ref.centroid_position[1] - ref.bounding_box[2]
    descr[0, 2] = 1.0 / (ref.intensity * counter)
    peak2d = subtrac_bkg_2d(data2d, background2d)
    sumation = add_2d(descr, peak2d, sumation)

  return sumation, thold

def fit_profile_2d(reflections, arr_proff, row, col, xmax, ymax):
  import math
  local_average = arr_proff[row][col][0]
  thold = arr_proff[row][col][1]

  len_tabl = len(arr_proff)
  descr = flex.double(flex.grid(1, 3))
  x_cuad_size = float(xmax) / len_tabl
  x_half_cuad_size = (x_cuad_size) / 2.0
  y_cuad_size = float(ymax) / len_tabl
  y_half_cuad_size = (y_cuad_size) / 2.0

  if_you_want_to_see_how_the_profiles_look = '''
  from matplotlib import pyplot as plt
  data2d = local_average.as_numpy_array()
  plt.imshow(data2d, interpolation = "nearest", cmap = plt.gray())
  plt.show()
  #'''

  for ref in reflections:
    if ref.is_valid() and ref.intensity < thold:


      x, y = ref.image_coord_px            # consider replasing with centroid pos

      if (x > x_half_cuad_size        and y > y_half_cuad_size and
          x < xmax - x_half_cuad_size and y < ymax - y_half_cuad_size):

        x_centr_of_cuad = col * x_cuad_size + x_half_cuad_size
        y_centr_of_cuad = row * y_cuad_size + y_half_cuad_size

        if x < x_centr_of_cuad and y < y_centr_of_cuad:
          tp_lf_pos = row - 1, col - 1
          tp_rg_pos =row - 1, col
          bt_lf_pos = row, col - 1
          bt_rg_pos = row, col
        elif x > x_centr_of_cuad and y < y_centr_of_cuad:
          tp_lf_pos = row - 1, col
          tp_rg_pos = row - 1, col + 1
          bt_lf_pos = row, col
          bt_rg_pos = row, col + 1
        elif x < x_centr_of_cuad and y > y_centr_of_cuad:
          tp_lf_pos = row, col - 1
          tp_rg_pos = row, col
          bt_lf_pos = row + 1, col - 1
          bt_rg_pos = row + 1, col
        else:
          tp_lf_pos = row, col
          tp_rg_pos = row, col + 1
          bt_lf_pos = row + 1, col
          bt_rg_pos = row + 1, col + 1

        tp_lf_average = arr_proff[tp_lf_pos[0]][tp_lf_pos[1]][0]
        tp_rg_average = arr_proff[tp_rg_pos[0]][tp_rg_pos[1]][0]
        bt_lf_average = arr_proff[bt_lf_pos[0]][bt_lf_pos[1]][0]
        bt_rg_average = arr_proff[bt_rg_pos[0]][bt_rg_pos[1]][0]



        tp_lf_x = tp_lf_pos[1] * x_cuad_size + x_half_cuad_size
        tp_lf_y = tp_lf_pos[0] * y_cuad_size + y_half_cuad_size
        dx = abs(tp_lf_x - x)
        dy = abs(tp_lf_y - y)
        tp_lf_dist = math.sqrt(dx * dx + dy * dy)

        tp_rg_x = tp_rg_pos[1] * x_cuad_size + x_half_cuad_size
        tp_rg_y = tp_rg_pos[0] * y_cuad_size + y_half_cuad_size
        dx = abs(tp_rg_x - x)
        dy = abs(tp_rg_y - y)
        tp_rg_dist = math.sqrt(dx * dx + dy * dy)

        bt_lf_x = bt_lf_pos[1] * x_cuad_size + x_half_cuad_size
        bt_lf_y = bt_lf_pos[0] * y_cuad_size + y_half_cuad_size
        dx = abs(bt_lf_x - x)
        dy = abs(bt_lf_y - y)
        bt_lf_dist = math.sqrt(dx * dx + dy * dy)

        bt_rg_x = bt_rg_pos[1] * x_cuad_size + x_half_cuad_size
        bt_rg_y = bt_rg_pos[0] * y_cuad_size + y_half_cuad_size
        dx = abs(bt_rg_x - x)
        dy = abs(bt_rg_y - y)
        bt_rg_dist = math.sqrt(dx * dx + dy * dy)

        max_dist = math.sqrt(
        (x_cuad_size * x_cuad_size) + (y_cuad_size * y_cuad_size)
        )
        tp_lf_contr = (max_dist - tp_lf_dist) / max_dist
        tp_rg_contr = (max_dist - tp_rg_dist) / max_dist
        bt_lf_contr = (max_dist - bt_lf_dist) / max_dist
        bt_rg_contr = (max_dist - bt_rg_dist) / max_dist

        total_contr = tp_lf_contr + tp_rg_contr + bt_lf_contr + bt_rg_contr

        re_scale = 1.0 / total_contr
        tp_lf_contr = tp_lf_contr * re_scale
        tp_rg_contr = tp_rg_contr * re_scale
        bt_lf_contr = bt_lf_contr * re_scale
        bt_rg_contr = bt_rg_contr * re_scale

        debugging_code = '''
        print "max_dist =", max_dist
        print "tp_lf_dist =", tp_lf_dist
        print "tp_rg_dist =", tp_rg_dist
        print "bt_lf_dist =", bt_lf_dist
        print "bt_rg_dist =", bt_rg_dist
        total_contr = tp_lf_contr + tp_rg_contr + bt_lf_contr + bt_rg_contr
        print "total_contr =", total_contr
        '''

        big_nrow = tp_lf_average.all()[0]
        if tp_rg_average.all()[0] > big_nrow:
          big_nrow = tp_rg_average.all()[0]
        if bt_lf_average.all()[0] > big_nrow:
          big_nrow = bt_lf_average.all()[0]
        if bt_rg_average.all()[0] > big_nrow:
          big_nrow = bt_rg_average.all()[0]

        big_ncol = tp_lf_average.all()[1]
        if tp_rg_average.all()[1] > big_ncol:
          big_ncol = tp_rg_average.all()[1]
        if bt_lf_average.all()[1] > big_ncol:
          big_ncol = bt_lf_average.all()[1]
        if bt_rg_average.all()[1] > big_ncol:
          big_ncol = bt_rg_average.all()[1]

        average = flex.double(flex.grid(big_nrow, big_ncol), 0)

        descr[0, 0] = float(tp_lf_average.all()[1])/2.0
        descr[0, 1] = float(tp_lf_average.all()[0])/2.0
        descr[0, 2] = float(tp_lf_contr)
        average = add_2d(descr, tp_lf_average, average)
        descr[0, 0] = float(tp_rg_average.all()[1])/2.0
        descr[0, 1] = float(tp_rg_average.all()[0])/2.0
        descr[0, 2] = float(tp_rg_contr)
        average = add_2d(descr, tp_rg_average, average)
        descr[0, 0] = float(bt_lf_average.all()[1])/2.0
        descr[0, 1] = float(bt_lf_average.all()[0])/2.0
        descr[0, 2] = float(bt_lf_contr)
        average = add_2d(descr, bt_lf_average, average)
        descr[0, 0] = float(bt_rg_average.all()[1])/2.0
        descr[0, 1] = float(bt_rg_average.all()[0])/2.0
        descr[0, 2] = float(bt_rg_contr)
        average = add_2d(descr, bt_rg_average, average)

        if_you_want_to_see_interpolated_profiles = '''
        if tp_lf_contr > 0.62 or tp_rg_contr > 0.62 \
        or bt_lf_contr > 0.62 or bt_rg_contr > 0.62:
          #from matplotlib import pyplot# as plt
          print "tp_lf_contr =", tp_lf_contr
          data2d = tp_lf_average.as_numpy_array()
          #pyplot.imshow(data2d, interpolation = "nearest", cmap = pyplot.gray())
          #pyplot.show()
          print "tp_rg_contr =", tp_rg_contr
          data2d = tp_rg_average.as_numpy_array()
          #pyplot.imshow(data2d, interpolation = "nearest", cmap = pyplot.gray())
          #pyplot.show()
          print "bt_lf_contr =", bt_lf_contr
          data2d = bt_lf_average.as_numpy_array()
          #pyplot.imshow(data2d, interpolation = "nearest", cmap = pyplot.gray())
          #pyplot.show()
          print "bt_rg_contr =", bt_rg_contr
          data2d = bt_rg_average.as_numpy_array()
          #pyplot.imshow(data2d, interpolation = "nearest", cmap = pyplot.gray())
          #pyplot.show()
          print "final averaged profile"
          data2d = average.as_numpy_array()
          data2d = average.as_numpy_array()
          #pyplot.imshow(data2d, interpolation = "nearest", cmap = pyplot.gray())
          #pyplot.show()
          print "____________________________"
        #'''


        if_you_want_to_see_how_the_profiles_look = '''
        from matplotlib import pyplot as plt
        data2d = average.as_numpy_array()
        plt.imshow(data2d, interpolation = "nearest", cmap = plt.gray())
        plt.show()
        #'''
      else:
        #print "in else"
        average = local_average



      shoebox = ref.shoebox
      mask = ref.shoebox_mask                               # may be needed soon
      background = ref.shoebox_background
      tmp_i = ref.intensity
      tmp_v = ref.intensity_variance
      ref.intensity = 0.0
      ref.intensity_variance = 0.0

      for i in range(shoebox.all()[0]):
        data2d = shoebox[i:i + 1, :, :]
        mask2d = mask[i:i + 1, :, :]                        # may be needed soon
        background2d = background[i:i + 1, :, :]
        try:
          data2d.reshape(flex.grid(shoebox.all()[1:]))
          mask2d.reshape(flex.grid(shoebox.all()[1:]))      # may be needed soon
          background2d.reshape(flex.grid(shoebox.all()[1:]))

        except:
          print "error reshaping flex-array"
          print "ref.bounding_box", ref.bounding_box
          break

        descr[0, 0] = ref.centroid_position[0] - ref.bounding_box[0]
        descr[0, 1] = ref.centroid_position[1] - ref.bounding_box[2]
        descr[0, 2] = 1.0 #/ (ref.intensity * counter)
        #fully_record = 'yes'
        if(ref.status == 0):
        #if(fully_record == 'yes'):
          vec_data = (tmp_i, tmp_v)
          a_mat_flx = flex.double(flex.grid(4, 4))
          b_vec_flx = flex.double(flex.grid(4, 1))
          ok_lg = fitting_2d_multile_var_build_mat(descr, data2d, background2d, \
                                        average, vec_data, a_mat_flx, b_vec_flx)

          #if ok_lg == 0:
          a_mat = a_mat_flx.as_scitbx_matrix()
          b_mat = b_vec_flx.as_scitbx_matrix()
          try:
            x_mat = a_mat.inverse() * b_mat
            k_abc_vec = x_mat.as_flex_double_matrix()
          except:
            print "fail to do profile fitting  <<<<<<<<"
            k_abc_vec=(0,0,0,0)
          #else:
          #  print "ok_lg != 0"
          #  k_abc_vec=(0,0,0,0)

          ref.intensity += k_abc_vec[0]
          ref.intensity_variance += k_abc_vec[1]
        else:
          I_R = fitting_2d_partials(descr, data2d, background2d, average, tmp_i)
          ref.intensity += I_R[0]
          ref.intensity_variance += I_R[1]



        var = sigma_2d(ref.intensity, mask2d, background2d)
        #reslt = sigma_2d(ref.intensity, mask2d, background2d)
        #ref.intensity += reslt[0]
        ref.intensity_variance += var



  return reflections
