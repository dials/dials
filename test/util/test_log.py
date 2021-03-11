import logging

import dials.util.log


def test_LoggingContext():

    # configure logging
    dials.util.log.config(verbosity=2)

    # get some loggers
    idx_logger = logging.getLogger("dials.algorithms.indexing")
    dials_logger = logging.getLogger("dials")

    # check the logging level is as expected
    assert idx_logger.getEffectiveLevel() == logging.DEBUG
    assert dials_logger.getEffectiveLevel() == logging.DEBUG

    # selectively change logging level and check this has worked
    with dials.util.log.LoggingContext(idx_logger, logging.ERROR):
        assert idx_logger.getEffectiveLevel() == logging.ERROR
        assert dials_logger.getEffectiveLevel() == logging.DEBUG

    # check logging levels are as they were before
    assert idx_logger.getEffectiveLevel() == logging.DEBUG
    assert dials_logger.getEffectiveLevel() == logging.DEBUG

    # now check we can pass logger as a string
    with dials.util.log.LoggingContext("dials.algorithms.indexing", logging.WARNING):
        assert idx_logger.getEffectiveLevel() == logging.WARNING
        assert dials_logger.getEffectiveLevel() == logging.DEBUG

    # check logging levels are as they were before
    assert idx_logger.getEffectiveLevel() == logging.DEBUG
    assert dials_logger.getEffectiveLevel() == logging.DEBUG
