from __future__ import annotations

import contextlib
import dataclasses
import faulthandler
import functools
import io
import os
import pathlib
import signal
import sys
from typing import Optional

import tabulate as _tabulate
import tqdm

from libtbx.utils import Sorry

# Coloured exit code support on windows
try:
    import colorama
except ImportError:
    pass
else:
    colorama.init()

__all__ = [
    "debug_console",
    "debug_context_manager",
    "progress",
    "show_mail_handle_errors",
    "Sorry",
    "tabulate",
]

# Define the default tablefmt in dials
tabulate = functools.partial(_tabulate.tabulate, tablefmt="psql")
functools.update_wrapper(tabulate, _tabulate.tabulate)

# Customisable progressbar decorator for iterators.
#
# Utilizes the progress bar from the tqdm package, with modified defaults:
#   - By default, resize when terminal is resized (dynamic-columns)
#   - By default, disables the progress bar for non-tty output
#   - By default, the progress bar will be removed after completion
#
# Usage:
#   >>> from dials.util import progress
#   >>> for i in progress(range(10)):
#   ...     ...
#
# See https://github.com/tqdm/tqdm for more in-depth usage and options.
progress = functools.partial(tqdm.tqdm, disable=None, dynamic_ncols=True, leave=False)
functools.update_wrapper(
    functools.partial(tqdm.tqdm, disable=None, dynamic_ncols=True, leave=False),
    tqdm.tqdm,
)


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

        embed(user_ns=namespace, colors="neutral")
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

    import multiprocessing
    import threading
    from datetime import datetime

    class DCM:
        def __init__(self, name, log_func):
            self._ocm = original_context_manager
            if name:
                self._name = f"{name}-{hash(original_context_manager)}"
            else:
                self._name = str(hash(original_context_manager))
            self._log_func = log_func

        def format_event(self, event):
            return f"{datetime.now().strftime('%H:%M:%S.%f')} {self._name}: {event}\n"

        def log(self, event):
            self._log_func(self.format_event(event))

        def __enter__(self, *args, **kwargs):
            try:
                raise Exception()
            except Exception:
                parent = sys.exc_info()[2].tb_frame.f_back
            call_code = f"{parent.f_code.co_filename}:{parent.f_lineno}"
            call_process = multiprocessing.current_process().name
            call_thread = threading.currentThread().getName()
            self.log(f"Knock {call_process}:{call_thread} {call_code}")
            z = self._ocm.__enter__(*args, **kwargs)
            self.log(f"Enter {call_process}:{call_thread} {call_code}")
            return z

        def __exit__(self, *args, **kwargs):
            call_process = multiprocessing.current_process().name
            call_thread = threading.currentThread().getName()
            self.log(f"Exit {call_process}:{call_thread}")
            self._ocm.__exit__(*args, **kwargs)
            self.log(f"Left {call_process}:{call_thread}")

    return DCM(name, log_func)


@contextlib.contextmanager
def enable_faulthandler():
    # Don't clobber someone elses faulthandler settings
    if not faulthandler.is_enabled():
        # SIGUSR2 not available on windows
        # The attached STDERR might not support what faulthandler wants
        with contextlib.suppress(AttributeError, io.UnsupportedOperation):
            faulthandler.enable()
            faulthandler.register(signal.SIGUSR2, all_threads=True)
    yield


@contextlib.contextmanager
def show_mail_on_error():
    try:
        yield
    except Exception as e:
        text = "Please report this error to dials-support@lists.sourceforge.net:"
        if issubclass(e.__class__, Sorry):
            raise

        if len(e.args) == 1:
            e.args = (f"{text} {e.args[0]}",)
        else:
            e.args = (text, *e.args)

        raise


@contextlib.contextmanager
def make_sys_exit_red():
    """Make sys.exit call messages red, if appropriate"""
    try:
        yield
    except SystemExit as e:
        # Don't colour if requested not to, or user formatted already
        if (
            isinstance(e.code, str)
            and "NO_COLOR" not in os.environ
            and sys.stderr.isatty()
            and "\033" not in e.code
        ):
            e.code = f"\033[31m{e.code}\033[0m"
        raise


@contextlib.contextmanager
def show_mail_handle_errors():
    """
    Handle showing errors to the user, including contact email, and colors.

    If a stack trace occurs, the user will be give the contact email to report
    the error to. If a sys.exit is caught, it will be converted to red, if
    appropriate.

    On Windows this additionally changes the stdout/stderr encoding to UTF8,
    which will only have an effect when the program output is redirected to
    a file or another process. This change will outlive the calling context.
    """
    if os.name == "nt":
        import dxtbx.util

        dxtbx.util.encode_output_as_utf8()
    with enable_faulthandler(), make_sys_exit_red(), show_mail_on_error():
        yield


@dataclasses.dataclass
class HTCondorJobClassAd:
    cpus_provisioned: Optional[int] = None
    memory_provisioned: Optional[float] = None


def parse_htcondor_job_classad(filepath: pathlib.Path) -> HTCondorJobClassAd:
    with filepath.open() as fh:
        memory_provisioned = None
        cpus_provisioned = None
        for line in fh.readlines():
            line = line.lower().strip()
            if line.startswith("memoryprovisioned") and "=" in line:
                memory_provisioned = float(line.split("=")[1])
            elif line.startswith("cpusprovisioned") and "=" in line:
                cpus_provisioned = int(line.split("=")[1])
        return HTCondorJobClassAd(
            cpus_provisioned=cpus_provisioned,
            memory_provisioned=memory_provisioned,
        )
