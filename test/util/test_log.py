import logging
from itertools import chain
from typing import Any, List

import dials.util.log
from dials.util.mp import batch_multi_node_parallel_map, multi_node_parallel_map


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


def test_easy_mp_logging():
    """
    Check the filtering of log records from easy_mp child processes.

    Child processes spawned/forked with libtbx.easy_mp use a logger that caches the
    raised messages with dials.util.log.CacheHandler, rather than emitting them.
    Here we check that the 'dials' loggers in the child processes inherit the threshold
    logging level of the 'dials' logger in the parent process.
    """
    test_log_message = "Here's a test log message."

    def log_something(_: Any) -> List[logging.LogRecord]:
        """
        Create a dummy info log record.

        A little dummy function to pass to the dials.util.mp parallel map functions,
        which simply logs a single info message.

        Args:
            _:  Absorb the expected argument and do nothing with it.

        Returns:
            The log records created during the calling of this function.
        """
        dials.util.log.config_simple_cached()
        logger = logging.getLogger("dials")
        logger.info(test_log_message)
        (handler,) = logger.handlers
        return handler.records

    # Check that the child processes cache INFO log messages when the parent process
    # logger's threshold level is the default (DEBUG).
    dials.util.log.config(verbosity=2)

    results = multi_node_parallel_map(
        log_something, range(5), njobs=2, nproc=2, cluster_method="multiprocessing"
    )
    results_batch = batch_multi_node_parallel_map(
        log_something,
        range(5),
        njobs=2,
        nproc=2,
        cluster_method="multiprocessing",
        callback=lambda _: None,
    )

    assert results[0][0].msg == test_log_message
    assert results_batch[0][0][0].msg == test_log_message

    # Check that the child processes do not cache INFO log messages when the parent
    # process logger's threshold level is WARNING.
    logging.getLogger("dials").setLevel(logging.WARNING)

    results = multi_node_parallel_map(
        log_something, range(5), njobs=2, nproc=2, cluster_method="multiprocessing"
    )
    results_batch = batch_multi_node_parallel_map(
        log_something,
        range(5),
        njobs=2,
        nproc=2,
        cluster_method="multiprocessing",
        callback=lambda _: None,
    )

    assert not any(chain.from_iterable(results))
    assert not any(chain.from_iterable(chain.from_iterable(results_batch)))
