from __future__ import absolute_import, division, print_function

import logging
import pytest

from dials.util.log import config_simple_stdout, LoggingContext


def test_LoggingContext():
    with pytest.deprecated_call():
        # configure logging
        config_simple_stdout()

    # get some loggers
    idx_logger = logging.getLogger("dials.algorithms.indexing")
    dials_logger = logging.getLogger("dials")

    # check the logging level is as expected
    assert idx_logger.getEffectiveLevel() == logging.DEBUG
    assert dials_logger.getEffectiveLevel() == logging.DEBUG

    # selectively change logging level and check this has worked
    with LoggingContext(idx_logger, logging.ERROR):
        assert idx_logger.getEffectiveLevel() == logging.ERROR
        assert dials_logger.getEffectiveLevel() == logging.DEBUG

    # check logging levels are as they were before
    assert idx_logger.getEffectiveLevel() == logging.DEBUG
    assert dials_logger.getEffectiveLevel() == logging.DEBUG

    # now check we can pass logger as a string
    with LoggingContext("dials.algorithms.indexing", logging.WARNING):
        assert idx_logger.getEffectiveLevel() == logging.WARNING
        assert dials_logger.getEffectiveLevel() == logging.DEBUG

    # check logging levels are as they were before
    assert idx_logger.getEffectiveLevel() == logging.DEBUG
    assert dials_logger.getEffectiveLevel() == logging.DEBUG
