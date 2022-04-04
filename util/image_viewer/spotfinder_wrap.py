from __future__ import annotations

import sys
import threading

import wx

from .slip_viewer.frame import chooser_wrapper as _chooser_wrapper
from .viewer_tools import ZeroMQEvent

try:
    import platform
    import resource

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


class spot_wrapper:
    def __init__(self, params):
        self.params = params
        self.frame = None
        self._zmq_thread = None

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

        # If ZMQ required, setup the endpoint to dispatch messages to wx
        if self.params.zmq_endpoint is not None:
            self._setup_zmq_endpoint(self.params.zmq_endpoint)

        #   debug_memory_usage()
        app.MainLoop()

        # Shut down the ZMQ endpoint if we created one
        if self.params.zmq_endpoint is not None:
            self._shutdown_zmq_endpoint()

    def load_image(self, filename):
        from dials.util.spotfinder_frame import create_load_image_event

        if self.frame is not None:
            create_load_image_event(self.frame, filename)

    def _setup_zmq_endpoint(self, endpoint: str) -> None:
        """Create and bind a ZMQ endpoint"""
        try:
            import zmq
        except ImportError:
            sys.exit("Error: ZMQ Endpoint requested but pyzmq not installed; aborting")
        print("ZeroMQ binding requested at", endpoint)
        self.zmq_context = zmq.Context()
        self.zmq_socket = self.zmq_context.socket(zmq.PULL)
        self.zmq_socket.bind(endpoint)
        self._zmq_thread = ZMQEventListener(self.zmq_socket, self.frame)
        self._zmq_thread.start()

    def _shutdown_zmq_endpoint(self):
        """Close the ZMQ endpoint"""
        print("Shutting down ZMQ endpoint")
        # Wait for the thread to end
        if self._zmq_thread is not None:
            self._zmq_thread.stop = True
            self._zmq_thread.join()

        self.zmq_socket.close()


class ZMQEventListener(threading.Thread):
    """Thread to listen for ZMQ events and dispatch them to the GUI.

    Args:
        socket: The ZeroMQ socket object to listen on
        target: The WX container to forward messages to
    """

    def __init__(self, socket, target, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.stop = False
        self.socket = socket
        self.target = target

    def run(self):
        while not self.stop:
            event_count = self.socket.poll(timeout=0.1)
            for _ in range(event_count):
                message = self.socket.recv_json()
                print("Received ZMQ Message: ", message)
                self._handle_message(message)

    def _handle_message(self, message):
        # Just attach the message to an event and send it on
        evt = ZeroMQEvent(message=message)
        wx.PostEvent(self.target, evt)
