from __future__ import absolute_import, division, print_function

import logging.config
import os
import six
import sys
import time

try:
    from dlstbx.util.colorstreamhandler import ColorStreamHandler
except ImportError:
    ColorStreamHandler = None

# https://stackoverflow.com/questions/25194864/python-logging-time-since-start-of-program/25196134#25196134
class ElapsedFormatter:
    """A formatter for log files that prepends messages with the elapsed time
    since initialisation and prefixes warning messages with 'WARNING:'"""

    def __init__(self):
        self.start_time = time.time()
        self.elapsed_msg = ""

    def format(self, record):
        elapsed_seconds = record.created - self.start_time
        elapsed_msg = "{:6.1f}: ".format(elapsed_seconds)
        indent = len(elapsed_msg)
        msg = record.getMessage()

        if record.levelno == logging.WARNING:
            msg = "WARNING: " + msg

        msg = msg.replace("\n", "\n" + " " * indent)
        if elapsed_msg == self.elapsed_msg:
            return " " * indent + msg
        else:
            self.elapsed_msg = elapsed_msg
            return elapsed_msg + msg


class LevelPrefixFormatter(logging.Formatter):
    def format(self, record):

        format_orig = self._fmt
        if record.levelno == logging.WARNING:
            self._fmt = "WARNING: %(msg)s"
        result = logging.Formatter.format(self, record)
        self._fmt = format_orig

        return result


def config(verbosity=0, logfile=None):
    """
    Configure the logging.

    :param verbosity: Verbosity level of log output. Possible values:
                        * 0: Info log output to stdout/logfile
                        * 1: Info & debug log output to stdout/logfile
    :type verbosity: int
    :param logfile: Filename for log output.  If False, no log file is written.
    :type logfile: str
    """

    if os.getenv("COLOURLOG") and ColorStreamHandler:
        console = ColorStreamHandler(sys.stdout)
    else:
        console = logging.StreamHandler(sys.stdout)

    fmt = LevelPrefixFormatter()
    console.setFormatter(fmt)
    dials_logger = logging.getLogger("dials")
    dials_logger.addHandler(console)

    if verbosity:
        loglevel = logging.DEBUG
    else:
        loglevel = logging.INFO

    if logfile:
        fh = logging.FileHandler(filename=logfile, mode="w")
        fh.setLevel(loglevel)
        fh.setFormatter(ElapsedFormatter())
        dials_logger.addHandler(fh)

    dials_logger.setLevel(loglevel)
    #   logging.getLogger("dxtbx").setLevel(logging.DEBUG)
    console.setLevel(loglevel)

    print_banner(use_logging=True)


class CacheHandler(logging.Handler):
    """A simple class to store log messages."""

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


class LoggingContext(object):
    # https://docs.python.org/3/howto/logging-cookbook.html#using-a-context-manager-for-selective-logging
    def __init__(self, logger, level=None):
        self.logger = (
            logging.getLogger(logger)
            if isinstance(logger, six.string_types)
            else logger
        )
        self.level = level

    def __enter__(self):
        if self.level is not None:
            self.old_level = self.logger.level
            self.logger.setLevel(self.level)

    def __exit__(self, et, ev, tb):
        if self.level is not None:
            self.logger.setLevel(self.old_level)
        # implicit return of None => don't swallow exceptions
