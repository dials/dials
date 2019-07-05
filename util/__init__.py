from __future__ import absolute_import, division, print_function

import contextlib
import sys
import warnings

from ._progress import progress  # noqa: F401, exported symbol

from libtbx.utils import Sorry


def debug_console():
    """Start python console at the current code point."""

    # use exception trick to pick up the current frame
    try:
        raise RuntimeError()
    except RuntimeError:
        frame = sys.exc_info()[2].tb_frame.f_back

    # evaluate commands in current namespace
    namespace = frame.f_globals.copy()
    namespace.update(frame.f_locals)

    try:
        # Start IPython console if IPython is available.
        from IPython import embed

        embed(user_ns=namespace)
    except ImportError:
        # Otherwise use basic python console
        import code

        code.interact(banner="=" * 80, local=namespace)


def debug_context_manager(original_context_manager, name="", log_func=None):
    """A wrapper to help debugging arbitrary context managers. For example
    instead of
      lock = threading.RLock():
    use
      lock = debug_context_manager(threading.RLock())
    and all calls
      with lock:
        ...
    will produce (hopefully) helpful debug output.

    :param original_context_manager: Some context manager to be wrapped.
    :param name: A name for the context manager, that will be printed in the
                 debug output.
    :param log_func: optional log function. If not specified debug output
                     will be printed to sys.stderr.
    :return: A context manager wrapping the original context manager.
    """
    if not log_func:

        def log_func(output):
            sys.stderr.write(output)
            sys.stderr.flush()

    from datetime import datetime
    import threading
    import multiprocessing

    class DCM(object):
        def __init__(self, name, log_func):
            self._ocm = original_context_manager
            if name:
                self._name = "%s-%s" % (name, str(hash(original_context_manager)))
            else:
                self._name = str(hash(original_context_manager))
            self._log_func = log_func

        def format_event(self, event):
            return "%s %s: %s\n" % (
                datetime.now().strftime("%H:%M:%S.%f"),
                self._name,
                event,
            )

        def log(self, event):
            self._log_func(self.format_event(event))

        def __enter__(self, *args, **kwargs):
            try:
                raise Exception()
            except Exception:
                parent = sys.exc_info()[2].tb_frame.f_back
            call_code = "%s:%s" % (parent.f_code.co_filename, parent.f_lineno)
            call_process = multiprocessing.current_process().name
            call_thread = threading.currentThread().getName()
            self.log("Knock %s:%s %s" % (call_process, call_thread, call_code))
            z = self._ocm.__enter__(*args, **kwargs)
            self.log("Enter %s:%s %s" % (call_process, call_thread, call_code))
            return z

        def __exit__(self, *args, **kwargs):
            call_process = multiprocessing.current_process().name
            call_thread = threading.currentThread().getName()
            self.log("Exit %s:%s" % (call_process, call_thread))
            z = self._ocm.__exit__(*args, **kwargs)
            self.log("Left %s:%s" % (call_process, call_thread))

    return DCM(name, log_func)


def halraiser(e):
    """
    Code:
      try:
        run()
      except Exception as e:
        halraiser(e)

    can be rewritten as
      with dials.util.show_mail_on_error():
        run()

    with the benefits of less boilerplate and that
    halraiser() does not appear in the traceback
    """
    warnings.warn(
        "halraiser() is deprecated. See https://github.com/dials/dials/commit/5c0d86142b8292f1017f5b0080b86917623dd4c9",
        DeprecationWarning,
        stacklevel=2,
    )

    text = "Please report this error to dials-support@lists.sourceforge.net:"

    if len(e.args) == 0:
        e.args = (text,)
    elif issubclass(e.__class__, Sorry):
        raise
    elif len(e.args) == 1:
        e.args = (text + " " + str(e.args[0]),)
    else:
        e.args = (text,) + e.args
    raise


@contextlib.contextmanager
def show_mail_on_error():
    try:
        yield
    except Exception as e:
        text = "Please report this error to dials-support@lists.sourceforge.net:"
        if len(e.args) == 0:
            e.args = (text,)
        elif issubclass(e.__class__, Sorry):
            raise
        elif len(e.args) == 1:
            e.args = (text + " " + str(e.args[0]),)
        else:
            e.args = (text,) + e.args
        raise
