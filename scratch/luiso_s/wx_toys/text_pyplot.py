import matplotlib.pyplot as plt
#import matplotlib as mpl
import numpy as np

fig = plt.figure()

ax = fig.add_subplot(1,1,1)

#ax.xaxis.set_ticklabels(["foo" , "bar", "ouch"])

plt.imshow(np.zeros( (50, 50), 'float'), interpolation = "nearest", vmin = 0, vmax = 55)
'''
print "get_view_interval =", ax.yaxis.get_view_interval()
print "get_data_interval =", ax.yaxis.get_data_interval()
print "get_label =", ax.yaxis.get_label()
print "get_majorticklabels =", ax.yaxis.get_majorticklabels()
print "get_label_text() =", ax.yaxis.get_label_text()
print "get_major_formatter() =", ax.yaxis.get_major_formatter()
'''
print "get_majorticklocs() =", ax.yaxis.get_majorticklocs()
'''
print "get_minor_ticks(numticks=None) =", ax.yaxis.get_minor_ticks(numticks=None)
print "get_minorticklocs() =", ax.yaxis.get_minorticklocs()
print "get_offset_text() =", ax.yaxis.get_offset_text()
print "get_scale() =", ax.yaxis.get_scale()
print "get_ticklabels(minor=False) =", ax.yaxis.get_ticklabels(minor=False)

print "get_ticklocs(minor=False) =", ax.yaxis.get_ticklocs(minor=False)
print "get_view_interval() =", ax.yaxis.get_view_interval()
print "get_units() =", ax.yaxis.get_units()
'''
labl = ax.yaxis.get_majorticklocs()
for pos in range(len(labl)):
  labl[pos] = labl[pos] + 3213.432423
ax.yaxis.set_ticklabels(labl)
#ax.yaxis.set_ticklabels(["Y_foo" , "Y_bar", "Y_ouch"])
plt.show()

