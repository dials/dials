def tst_dxtbx():
    import libtbx.load_env
    import os
    dials_regression = libtbx.env.dist_path('dials_regression')
    from boost.python import streambuf
    from dxtbx import read_uint16
    from dxtbx.format.Registry import Registry

    for directory, image in [('SLS_X06SA', 'mar225_2_001.img')]:
        file_path = os.path.join(dials_regression, 'image_examples',
                                 directory, image)
        format = Registry.find(file_path)
        i = format(file_path)
        size = i.get_detector().get_image_size()

    print 'OK'

tst_dxtbx()
