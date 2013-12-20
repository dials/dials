#from dxtbx.model import Panel

#dl = (0.3798037784837018, 0.9158092591648166, 105.08859002092663,
#      -0.384891609331586, 0.2847248712008545, -118.4568375432652,
#      -0.8411941146102417, 0.2832158912989938, -192.95017781758185)

#dl1 = dl[0], dl[3], dl[6]
#dl2 = dl[1], dl[4], dl[7]
#dl0 = dl[2], dl[5], dl[8]

#d1 = (0.379803835664, -0.384891594186, -0.841194095175)
#d2 = (0.915801577474, 0.284722451667, 0.283213456118)
#d0 = (105.088371189, -118.45645418, -192.950453736)

#from scitbx import matrix
#d1 = matrix.col(d1)
#d2 = matrix.col(d2)

#print d1.dot(d2)
#p = Panel()
#p.set_local_frame(dl1, dl2, dl0)
#p.set_frame(d1, d2, d0)

#d11 = p.get_fast_axis()
#d22 = p.get_slow_axis()
#d00 = p.get_origin()

#print d11
#print d22
#print d00


from dxtbx.model import Panel

dp = (1.0, 0.0, -20.955, 0.0, -1.0, 10.615, 0.0, 0.0, 0.0)

dl = (-0.001825172553060027, -0.9999981947571188, -37.76001244640159, 0.9999973206373582, -0.0018259228624977202, 19.59161224494454, -0.0014238901839527113, -0.0005258214562373528, -150.0059)

dl1 = dl[0], dl[3], dl[6]
dl2 = dl[1], dl[4], dl[7]
dl0 = dl[2], dl[5], dl[8]

dp1 = dp[0], dp[3], dp[6]
dp2 = dp[1], dp[4], dp[7]
dp0 = dp[2], dp[5], dp[8]

d1 = (-0.00182517255306, 0.999997320637, -0.00142389018396)
d2 = (0.999998194757, 0.0018259228625, 0.000525821456245)
d0 = (-48.3367467929, -1.3827137802, -149.981643976)

p = Panel()
p.set_parent_frame(dp1, dp2, dp0)
p.set_local_frame(dl1, dl2, dl0)
p.set_frame(d1, d2, d0)


#print repr([10000000000*pp for pp in p.get_parent_d_matrix()])
#for pp in p.get_parent_d_matrix():
#  print '{0:.16f}'.format(pp)

#
#p.set_local_frame(dl1, dl2, dl0)

#d_matrix1 = p.get_d_matrix()

#p = Panel()
#p.set_frame(dl1, dl2, dl0)
#d_matrix2 = p.get_d_matrix()

#for d1, d2 in zip(d_matrix1, d_matrix2):
#  print d1 - d2
