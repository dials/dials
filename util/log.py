from __future__ import absolute_import, division, print_function

import logging.config
import os
import sys
import warnings

try:
    from dlstbx.util.colorstreamhandler import ColorStreamHandler
except ImportError:
    ColorStreamHandler = None


def config(verbosity=0, name=None, info=None, debug=None, logfile=None):
    """
    Configure the logging.

    :param verbosity: Verbosity level of log output. Possible values:
                        * 0: Info log output to stdout/logfile
                        * 1: Info & debug log output to stdout/logfile
    :type verbosity: int
    :param logfile: Filename for log output.  If False, no log file is written.
    :type logfile: str
    """

    if info:
        warnings.warn(
            "info= parameter is deprecated, use logfile=",
            DeprecationWarning,
            stacklevel=2,
        )
    if debug:
        warnings.warn(
            "debug= parameter is deprecated, use logfile= and verbosity=",
            DeprecationWarning,
            stacklevel=2,
        )
    if name:
        warnings.warn("name= parameter is deprecated", DeprecationWarning, stacklevel=2)

    if os.getenv("COLOURLOG") and ColorStreamHandler:
        console = ColorStreamHandler(sys.stdout)
    else:
        console = logging.StreamHandler(sys.stdout)

    dials_logger = logging.getLogger("dials")
    dials_logger.addHandler(console)

    if verbosity:
        loglevel = logging.DEBUG
    else:
        loglevel = logging.INFO

    logfilename = logfile or info or debug
    if logfilename:
        fh = logging.FileHandler(filename=logfilename, mode="w")
        fh.setLevel(loglevel)
        dials_logger.addHandler(fh)

    dials_logger.setLevel(loglevel)
    #   logging.getLogger("dxtbx").setLevel(logging.DEBUG)
    console.setLevel(loglevel)

    print_banner(use_logging=True)


def config_simple_stdout(name="dials"):
    warnings.warn(
        "config_simple_stdout is deprecated, use config",
        DeprecationWarning,
        stacklevel=2,
    )

    """
    Configure the logging to just go to stdout
    """

    # Configure the logging
    logging.config.dictConfig(
        {
            "version": 1,
            "disable_existing_loggers": False,
            "formatters": {"standard": {"format": "%(message)s"}},
            "handlers": {
                "stream": {
                    "level": "DEBUG",
                    "class": "logging.StreamHandler",
                    "formatter": "standard",
                    "stream": "ext://sys.stdout",
                }
            },
            "loggers": {
                name: {"handlers": ["stream"], "level": "DEBUG", "propagate": True}
            },
        }
    )

    print_banner(use_logging=True)


class CacheHandler(logging.Handler):
    """ A simple class to store log messages. """

    def __init__(self):
        """
        Initialise the handler
        """
        super(CacheHandler, self).__init__()
        self._messages = []

    def emit(self, record):
        """
        Emit the message to a list

        :param record: The log record
        """
        self._messages.append(record)

    def messages(self):
        return self._messages


def config_simple_cached():
    """
    Configure the logging to use a cache.
    """

    # Configure the logging
    logging.config.dictConfig(
        {
            "version": 1,
            "disable_existing_loggers": False,
            "handlers": {
                "cache": {"level": "DEBUG", "class": "dials.util.log.CacheHandler"}
            },
            "loggers": {
                "dials": {"handlers": ["cache"], "level": "DEBUG", "propagate": True}
            },
        }
    )


_banner = (
    "DIALS (2018) Acta Cryst. D74, 85-97. https://doi.org/10.1107/S2059798317017235"
)
_banner_printed = False


def print_banner(force=False, use_logging=False):
    global _banner_printed
    if _banner_printed and not force:
        return
    if os.getenv("DIALS_NOBANNER"):
        return
    _banner_printed = True

    if use_logging:
        logging.getLogger("dials").info(_banner)
    else:
        print(_banner)
