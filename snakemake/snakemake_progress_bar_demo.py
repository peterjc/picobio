#!/usr/bin/env python3
"""Demonstration of calling a snakemake workflow with a progress bar.

Written and tested using snakemake 8.20.6 under macOS.

Currently the snakemake API doesn't have any obvious way
to get callbacks or an iterator approach to running a
workflow which would allow direct updates to a progress
bar. Improvements to their logging system may allow this?

Instead, this demonstrate running snakemake on a subprocess,
and monitoring the creation of the expected output files
as a proxy to update a progress bar. This works, but would
put some additional load on the file system.

This uses the rich library's progress bar, but the same idea
would work with another library like tqdm. We must explicitly
update the progress bar whenever a new output file is found.
"""

from multiprocessing import Process
from pathlib import Path

from rich.progress import Progress  # or use tqdm, or ...

from snakemake.api import DAGSettings
from snakemake.api import ResourceSettings
from snakemake.api import SnakemakeApi
from snakemake.settings.types import OutputSettings
from snakemake.settings.types import Quietness

inputs = Path(".").glob("*.fna")
targets = [str(_) + ".md5" for _ in inputs]


def black_box(output_files):
    """Black-box function which generates known files (here snakemake)."""
    snakefile = Path("demo.smk")
    with SnakemakeApi(OutputSettings(quiet={Quietness.ALL})) as snakemake_api:
        workflow_api = snakemake_api.workflow(
            snakefile=snakefile,
            resource_settings=ResourceSettings(),
            # config_settings=ConfigSettings(config=config_args),
            # workdir=workdir,
        )
        dag_api = workflow_api.dag(
            dag_settings=DAGSettings(targets=output_files),
        )
        dag_api.unlock()
        dag_api.execute_workflow()


def with_progress_bar(function, output_files, interval=0.5):
    """Run given function via subprocess with a progress bar.

    The function must accept a single argument, the given file list.
    The appearance of those files on disk is used to update the progress
    bar. This runs the function in a process via multiprocessing, and
    returns the process exit code (should be zero for success).
    """
    pending = [Path(_) for _ in targets]
    p = Process(target=function, args=(output_files,))
    p.start()
    with Progress() as progress:
        task = progress.add_task("Snakemake...", total=len(pending))
        while pending:
            p.join(interval)
            for t in pending[:]:
                if t.is_file():
                    print(f"Done: {t}")
                    pending.remove(t)
                    progress.update(task, advance=1)
            if p.exitcode is not None:
                # Should be finished, but was it success or failure?
                pending = []  # to break the loop
    p.join()  # Should be immediate as should have finished
    assert not p.is_alive()
    print(f"Snakemake return code {p.exitcode}")
    return p.exitcode


if __name__ == "__main__":
    # black_box(targets)
    with_progress_bar(black_box, targets)
