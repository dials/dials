from __future__ import division

if __name__ == '__main__':

    from extract_reflection_profiles import load_cbf_image_volume
    from matplotlib import pylab
    from matplotlib import cm

    image_volume = load_cbf_image_volume('/home/upc86896/Projects/data/300k/ximg2700*.cbf')

    image_size = image_volume.shape
    image_size_ratio = image_size[2] / image_size[1]
    figure_size = (6, 6 / image_size_ratio)

    for i in range(0, 9):
        fig = pylab.figure(figsize=figure_size, dpi=300)
        plt = pylab.imshow(image_volume[i,:,:], vmin=0, vmax=1000, cmap=cm.Greys_r)
        plt.axes.get_xaxis().set_ticks([])
        plt.axes.get_yaxis().set_ticks([])
        fig.savefig('/home/upc86896/Documents/image_volume_{0}.tiff'.format(i),
            bbox_inches='tight', pad_inches=0)
