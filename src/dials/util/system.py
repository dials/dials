from __future__ import annotations

import logging
import math
import os
import sys

import psutil

__all__ = ("cpu_count", "memory_limit", "CPU_COUNT", "MEMORY_LIMIT")


logger = logging.getLogger(__name__)


def _try_extract_cgroup_cpu_quota():
    # cgroup v1
    # The directory name isn't standardized across linux distros, check both
    for dirname in ["cpuacct,cpu", "cpu,cpuacct"]:
        try:
            with open("/sys/fs/cgroup/%s/cpu.cfs_quota_us" % dirname) as f:
                quota = int(f.read())
            with open("/sys/fs/cgroup/%s/cpu.cfs_period_us" % dirname) as f:
                period = int(f.read())
            return quota, period
        except Exception:
            pass

    # cgroup v2
    try:
        with open("/proc/self/cgroup") as f:
            group_path = f.read().strip().split(":")[-1]
        if not group_path.endswith("/"):
            group_path = f"{group_path}/"
        with open("/sys/fs/cgroup%scpu.max" % group_path) as f:
            quota, period = map(int, f.read().split(" "))
            return quota, period
    except Exception:
        pass

    # No cgroup CPU quota found
    return None, None


def cpu_count() -> int:
    """Get the available CPU count for this system.

    Takes the minimum value from the following locations:

    - Total system cpus available on the host.
    - CPU Affinity (if set)
    - Cgroups limit (if set)

    If the number of CPUs cannot be determined by any of these means, return 1.
    """
    count: int | None = os.cpu_count()

    # Check CPU affinity if available
    try:
        if affinity := psutil.Process().cpu_affinity():
            affinity_count = len(affinity)
            if affinity_count > 0:
                count = min(count, affinity_count) if count else affinity_count
    except Exception:
        pass

    # Check cgroups if available
    if sys.platform == "linux":
        quota, period = _try_extract_cgroup_cpu_quota()
        if quota is not None and period is not None:
            # We round up on fractional CPUs
            cgroups_count = math.ceil(quota / period)
            if cgroups_count > 0:
                count = min(count, cgroups_count) if count else cgroups_count

    return count or 1


def memory_limit() -> int:
    """Get the memory limit (in bytes) for this system.

    Takes the minimum value from the following locations:

    - Available system host memory
    - Cgroups limit (if set)
    - RSS rlimit (if set)
    """
    limit = psutil.virtual_memory().available

    # Check cgroups if available
    # Note: can't use LINUX and WINDOWS constants as they upset mypy
    if sys.platform == "linux":
        path_used = None
        for path in [
            "/sys/fs/cgroup/memory/memory.limit_in_bytes",  # cgroups v1 hard limit
            "/sys/fs/cgroup/memory/memory.soft_limit_in_bytes",  # cgroups v1 soft limit
            "/sys/fs/cgroup/memory.max",  # cgroups v2 hard limit
            "/sys/fs/cgroup/memory.high",  # cgroups v2 soft limit
            "/sys/fs/cgroup/memory.low",  # cgroups v2 softest limit
        ]:
            try:
                with open(path) as f:
                    cgroups_limit = int(f.read())
                if cgroups_limit > 0:
                    path_used = path
                    limit = min(limit, cgroups_limit)
            except Exception:
                pass
        if path_used:
            logger.debug(
                "Setting system memory limit based on cgroup value defined in %s",
                path_used,
            )

    # Check rlimit if available
    if sys.platform != "win32":
        try:
            import resource

            hard_limit = resource.getrlimit(resource.RLIMIT_RSS)[1]
            if 0 < hard_limit < limit:
                logger.debug(
                    "Limiting system memory based on RLIMIT_RSS to %s", hard_limit
                )
                limit = hard_limit
        except (ImportError, OSError):
            pass

    return limit


CPU_COUNT = cpu_count()
MEMORY_LIMIT = memory_limit()
