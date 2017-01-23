from __future__ import absolute_import, division
import os

import libtbx.load_env

have_dials_regression = libtbx.env.has_module("dials_regression")
if have_dials_regression:
  dials_regression = libtbx.env.find_in_repositories(
    relative_path="dials_regression",
    test=os.path.isdir)


def exercise_spot_counts_per_image_plot():
  from libtbx import easy_pickle
  from dials.util import ascii_art
  from libtbx.test_utils import show_diff
  if not have_dials_regression:
    print 'Skipping exercise_spot_counts_per_image_plot(): dials_regression not available'
    return
  data_dir = os.path.join(dials_regression, "indexing_test_data", "i04_weak_data")
  pickle_path = os.path.join(data_dir, "full.pickle")

  refl = easy_pickle.load(pickle_path)
  output = ascii_art.spot_counts_per_image_plot(refl)

  expected_output = '''\
116082 spots found on 540 images (max 2072 / bin)
* *      *                     *    * * *  * **    *       *
************************************************************
************************************************************
************************************************************
************************************************************
************************************************************
************************************************************
************************************************************
************************************************************
************************************************************
1                         image                          540'''
  output = '\n'.join(line.rstrip() for line in output.splitlines())
  assert not show_diff(output, expected_output)


  output = ascii_art.spot_counts_per_image_plot(
    refl, char='o', width=80, height=15)

  expected_output = '''\
116082 spots found on 540 images (max 1556 / bin)
o o         o o                           o o   o  o  o  oo                   oo
ooooooooooo o o ooo o oooo ooooooo ooo oo ooo ooo ooo oooooooooo ooooooo ooo ooo
oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
1                                   image                                    540'''
  output = '\n'.join(line.rstrip() for line in output.splitlines())
  assert not show_diff(output, expected_output)

  output = ascii_art.spot_counts_per_image_plot(
    refl, char='#', width=7, height=10)

  expected_output = '''\
116082 spots found on 540 images (max 16736 / bin)
#######
#######
#######
#######
#######
#######
#######
#######
#######
#######
1   540'''
  output = '\n'.join(line.rstrip() for line in output.splitlines())
  assert not show_diff(output, expected_output)

  output = ascii_art.spot_counts_per_image_plot(
    refl.select(refl['xyzobs.px.value'].parts()[2] <= 9.5), char='#', width=10, height=15)
  expected_output = '''\
2065 spots found on 10 images (max 342 / bin)
#
#
#
#
#
# #
# # #  ###
# #### ###
##########
##########
##########
##########
##########
##########
##########
1 image 10'''
  output = '\n'.join(line.rstrip() for line in output.splitlines())
  assert not show_diff(output, expected_output)


def run():
  exercise_spot_counts_per_image_plot()
  print 'OK'


if __name__ == '__main__':
  run()
