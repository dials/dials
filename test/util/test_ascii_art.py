from __future__ import absolute_import, division, print_function

import os

def test_spot_counts_per_image_plot(dials_regression):
  from libtbx import easy_pickle
  from dials.util import ascii_art
  data_dir = os.path.join(dials_regression, "indexing_test_data", "i04_weak_data")
  pickle_path = os.path.join(data_dir, "full.pickle")

  refl = easy_pickle.load(pickle_path)
  output = ascii_art.spot_counts_per_image_plot(refl)

  expected_output = '''\
116082 spots found on 540 images (max 2075 / bin)
* *      *        *                 * * *  * **            *
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
  assert output == expected_output

  output = ascii_art.spot_counts_per_image_plot(
    refl, char='o', width=80, height=15)

  expected_output = '''\
116082 spots found on 540 images (max 1627 / bin)
                                                         o
o oo   o  o o o       o o o o o  o  o     o o   o  o oo  oo o o   o o   o o   oo
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
  assert output == expected_output

  output = ascii_art.spot_counts_per_image_plot(
    refl, char='#', width=7, height=10)

  expected_output = '''\
116082 spots found on 540 images (max 16752 / bin)
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
  assert output == expected_output

  output = ascii_art.spot_counts_per_image_plot(
    refl.select(refl['xyzobs.px.value'].parts()[2] <= 9.5), char='#', width=10, height=15)
  expected_output = '''\
2065 spots found on 10 images (max 361 / bin)
#
#
#
#
#
#        #
#   #    #
## ##   ##
##########
##########
##########
##########
##########
##########
##########
1 image 10'''
  output = '\n'.join(line.rstrip() for line in output.splitlines())
  assert output == expected_output
