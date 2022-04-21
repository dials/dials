# LIBTBX_SET_DISPATCHER_NAME cluster.dials.exec

from __future__ import annotations

import pickle

import dials.util


def get_cwd():
    """
    Get the current working directory
    """
    import sys

    return sys.argv[1]


def get_tid():
    """
    Get the task id
    """
    import os

    # FIXME At the moment, there is no portable way to know the task id through
    # drmaa. This is really annoying. So we have to use the SGE_TASK_ID.
    # Therefore, this will only work for SGE. We can probably add support for
    # other systems as and when needed.
    if "SGE_TASK_ID" in os.environ:
        return os.environ["SGE_TASK_ID"]
    else:
        raise KeyError("Could not find task id")


@dials.util.show_mail_handle_errors()
def run(_=None):
    import traceback
    from os.path import exists, join
    from time import sleep

    # Get the task id and the current working directory
    tid = get_tid()
    cwd = get_cwd()

    # Set the paths
    input_fn = join(cwd, f"{tid}.input")
    output_fn = join(cwd, f"{tid}.output")

    # Wait until it exists
    while not exists(input_fn):
        sleep(1)

    # Try to run the function, otherwise return an exception
    try:
        with open(input_fn, "rb") as infile:
            function, element = pickle.load(infile)

        result = function(element)
    except Exception as e:
        e.args = [traceback.format_exc()]
        result = e

    # Dump the result
    with open(output_fn, "wb") as outfile:
        pickle.dump(result, outfile, protocol=pickle.HIGHEST_PROTOCOL)


if __name__ == "__main__":
    run()
