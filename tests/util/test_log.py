from __future__ import annotations

import logging
from typing import Any, List

import pytest

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
    test_log_message = "Here's a test log message."
    dials.util.log.config_simple_cached()
    logger = logging.getLogger("dials")
    logger.info(test_log_message)
    (handler,) = logger.handlers
    return handler.records


@pytest.mark.xfail(
    "os.name == 'nt'",
    reason="https://github.com/dials/dials/issues/1650",
    raises=AttributeError,
)
def test_cached_log_records(caplog):

    # Generate some cached log messages in easy_mp child processes.
    results = multi_node_parallel_map(
        log_something,
        iterable=range(5),
        njobs=2,
        nproc=2,
        cluster_method="multiprocessing",
    )
    # Get all the resulting log records in a single flattened list.
    results = [record for records in results for record in records]

    results_batch = batch_multi_node_parallel_map(
        log_something,
        range(5),
        njobs=2,
        nproc=2,
        cluster_method="multiprocessing",
        callback=lambda _: None,
    )
    # Get all the resulting log records in a single flattened list.
    results_batch = [
        record for batch in results_batch for records in batch for record in records
    ]

    # Check that re-handling the messages in a logger in this process with a
    # threshold severity of WARNING results in no log records being emitted.
    with dials.util.log.LoggingContext("dials", logging.WARNING):
        dials.util.log.rehandle_cached_records(results)
        assert not caplog.records

        dials.util.log.rehandle_cached_records(results_batch)
        assert not caplog.records

    # Check that re-handling the messages in a logger in this process with a
    # threshold severity of INFO results in all the log records being emitted.
    with dials.util.log.LoggingContext("dials", logging.INFO):
        dials.util.log.rehandle_cached_records(results)
        assert caplog.records == results
        caplog.clear()

        dials.util.log.rehandle_cached_records(results_batch)
        assert caplog.records == results_batch
        caplog.clear()
