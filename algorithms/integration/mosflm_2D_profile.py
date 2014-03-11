from __future__ import division
from dials.model.data import Reflection, ReflectionList
from dials.algorithms.integration import add_2d, subtrac_bkg_2d, \
              fitting_2d_partials, fitting_2d_multile_var_build_mat, sigma_2d

from dials.array_family import flex

def make_2d_profile(reflection_pointers, ref_table_in):
  '''
  t_shoebox = flex.shoebox(num_ref)
  t_intensity = flex.double(num_ref)
  t_xyzobs = flex.vec3_double(num_ref)
  '''
  #t_xyzcal = flex.vec3_double(num_ref)

  t_intensity = ref_table_in['intensity.raw.value']
  max_i_01 = 0.0
  for ref in reflection_pointers:
    #if ref.is_valid():
      if t_intensity[ref] > max_i_01:
        max_i_01 = t_intensity[ref]
  print "max_i_01 =", max_i_01

  max_i = 0.0
  for ref in reflection_pointers:
    #if ref.is_valid():
      if t_intensity[ref] > max_i and t_intensity[ref] < max_i_01 * 0.95:
        max_i = t_intensity[ref]
  thold = 0.5 * max_i
  print "max_i =", max_i
  print "thold =", thold


  select_pointers = []
  for ref in reflection_pointers:
    #if ref.is_valid():
      if t_intensity[ref] > thold and t_intensity[ref] < max_i:
        select_pointers.append(ref)
  counter = 0
  print "len(select_pointers) =", len(select_pointers)
  big_nrow = 0
  big_ncol = 0
  t_shoebox = ref_table_in['shoebox']
  for ref in select_pointers:
    local_nrow = t_shoebox[ref].data.all()[1]
    local_ncol = t_shoebox[ref].data.all()[2]
    if local_nrow > big_nrow:
      big_nrow = local_nrow
    if local_ncol > big_ncol:
      big_ncol = local_ncol
    counter += 1
  print big_nrow, big_ncol
  big_nrow = big_nrow * 2 + 1
  big_ncol = big_ncol * 2 + 1

  t_xyzobs = ref_table_in['xyzobs.px.value']
  t_bbox = ref_table_in['bbox']

  sumation = flex.double(flex.grid(big_nrow, big_ncol), 0)
  descr = flex.double(flex.grid(1, 3), 0)
  for ref in select_pointers:


    shoebox = t_shoebox[ref].data
    background = t_shoebox[ref].background

    data2d = shoebox[0:1, :, :]
    background2d = background[0:1, :, :]

    data2d.reshape(flex.grid(shoebox.all()[1:]))
    background2d.reshape(flex.grid(background.all()[1:]))

    #mask =t_shoebox[ref].mask                      # may be needed soon
    #mask2d = mask[0:1, :, :]                       # may be needed soon
    #mask2d.reshape(flex.grid(mask.all()[1:]))      # may be needed soon

    cntr_pos = t_xyzobs[ref]
    bnd_box = t_bbox[ref]

    descr[0, 0] = cntr_pos[0] - bnd_box[0]
    descr[0, 1] = cntr_pos[1] - bnd_box[2]
    descr[0, 2] = 1.0 / (t_intensity[ref] * counter)
    peak2d = subtrac_bkg_2d(data2d, background2d)
    sumation = add_2d(descr, peak2d, sumation)



  #if_you_want_to_see_how_the_profiles_look = '''
  from matplotlib import pyplot as plt
  data2d_np = sumation.as_numpy_array()
  plt.imshow(data2d_np, interpolation = "nearest", cmap = plt.gray())
  plt.show()
  #'''

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
          tmp_scale = tmp_i
          a_mat_flx = flex.double(flex.grid(4, 4))
          b_vec_flx = flex.double(flex.grid(4, 1))
          ok_lg = fitting_2d_multile_var_build_mat(descr, data2d, background2d, \
                                        average, tmp_scale, a_mat_flx, b_vec_flx)

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
