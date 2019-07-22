from __future__ import absolute_import, division, print_function

from .slip_viewer.frame import chooser_wrapper as _chooser_wrapper

try:
    import resource
    import platform

    def debug_memory_usage():
        # getrusage returns kb on linux, bytes on mac
        units_per_mb = 1024
        if platform.system() == "Darwin":
            units_per_mb = 1024 * 1024
        print(
            "Memory usage: %.1f MB"
            % (int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss) / units_per_mb)
        )


except ImportError:

    def debug_memory_usage():
        pass


class chooser_wrapper(_chooser_wrapper):
    def show_header(self):
        #   debug_memory_usage()
        pass

    def __eq__(self, other):
        return (
            hasattr(other, "image_set")
            and self.image_set is other.image_set
            and hasattr(other, "index")
            and self.index == other.index
        )

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash((id(self.image_set), self.index))


class spot_wrapper(object):
    def __init__(self, params):
        self.params = params
        self.frame = None

    def display(self, experiments, reflections):
        import wx
        from dials.util.image_viewer.spotfinder_frame import SpotFrame

        app = wx.App()

        self.frame = SpotFrame(
            None,
            -1,
            "X-ray image display",
            size=(800, 720),
            pos=(100, 100),
            params=self.params,
            experiments=experiments,
            reflections=reflections,
        )
        self.frame.SetSize((1024, 780))
        self.frame.Show()

        imagesets = self.frame.imagesets
        for imageset in imagesets:
            for idx in range(len(imageset.indices())):
                self.frame.add_file_name_or_data(chooser_wrapper(imageset, idx))
        # Make sure we load the first image as the current one
        self.frame.load_image(chooser_wrapper(imagesets[0], 0))
        #   debug_memory_usage()
        app.MainLoop()

    def load_image(self, filename):
        from dials.util.spotfinder_frame import create_load_image_event

        if self.frame is not None:
            create_load_image_event(self.frame, filename)
