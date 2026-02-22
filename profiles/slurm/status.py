#!/usr/bin/env python3
"""
status.py — SLURM job status checker for Snakemake --cluster-status

Snakemake calls this script with a job ID to determine whether a SLURM job
is still running, succeeded, or failed. Without this, Snakemake relies solely
on checking whether output files appeared, which is unreliable on networked
filesystems (NFS/Lustre metadata lag).

Usage by Snakemake (automatic):
    python profiles/slurm/status.py <jobid>

Returns one of: "running", "success", "failed"
"""

import subprocess
import sys
import time


# States that mean the job is still in progress
RUNNING_STATES = {
    "PENDING", "RUNNING", "REQUEUED", "SUSPENDED",
    "CONFIGURING", "COMPLETING", "RESIZING",
}

# States that mean the job succeeded
SUCCESS_STATES = {
    "COMPLETED",
}

# Everything else is a failure (FAILED, CANCELLED, TIMEOUT, OUT_OF_MEMORY,
# NODE_FAIL, PREEMPTED, BOOT_FAIL, DEADLINE, etc.)


def get_job_state(jobid: str, retries: int = 3) -> str:
    """Query sacct for the job state, with retries for transient failures."""
    for attempt in range(retries):
        try:
            result = subprocess.run(
                ["sacct", "-j", jobid, "--format=State", "--noheader", "--parsable2"],
                capture_output=True,
                text=True,
                timeout=30,
            )
            if result.returncode == 0 and result.stdout.strip():
                # sacct can return multiple lines (job + job steps).
                # Take the first non-empty state from the main job line.
                for line in result.stdout.strip().split("\n"):
                    state = line.strip().split("|")[0].strip()
                    if state:
                        return state
        except (subprocess.TimeoutExpired, FileNotFoundError):
            pass

        # Brief pause before retry
        if attempt < retries - 1:
            time.sleep(2 * (attempt + 1))

    # If sacct keeps failing, try squeue as fallback
    try:
        result = subprocess.run(
            ["squeue", "-j", jobid, "--format=%T", "--noheader"],
            capture_output=True,
            text=True,
            timeout=15,
        )
        if result.returncode == 0 and result.stdout.strip():
            return result.stdout.strip().split("\n")[0].strip()
    except (subprocess.TimeoutExpired, FileNotFoundError):
        pass

    # If we truly can't determine state, report failed to be safe
    # (Snakemake will rerun the job)
    return "UNKNOWN"


def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <jobid>", file=sys.stderr)
        sys.exit(1)

    jobid = sys.argv[1].strip()
    state = get_job_state(jobid)

    if state in RUNNING_STATES:
        print("running")
    elif state in SUCCESS_STATES:
        print("success")
    else:
        # Log the actual state for debugging
        print(f"failed (SLURM state: {state})", file=sys.stderr)
        print("failed")


if __name__ == "__main__":
    main()
