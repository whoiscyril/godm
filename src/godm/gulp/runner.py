from pathlib import Path
from ..models import Settings, Structure
from ..config import read_inputs
from ..gulp.inputs import read_library, write_master_input, write_ml_reference_input
from ..gulp.outputs import read_output, read_restart
import subprocess
from concurrent.futures import ThreadPoolExecutor


def run_gulp_batch(
    input_paths: list[Path],
    settings: Settings,
    ncpus: int,
    num_workers: int,
    launcher: str,
) -> list[Path]:
    """Perform concurrent GULP runs"""
    if ncpus < 1:
        raise ValueError("ncpu must be positive")
    if num_workers < 1:
        raise ValueError("num_workers must be positive")

    with ThreadPoolExecutor(max_workers=num_workers) as executor:
        futures = [
            executor.submit(
                run_gulp,
                ncpus=ncpus,
                input_path=path,
                directory=path.parent,
                executable=settings.gulp_path,
                launcher=launcher,
                timeout=None,
            )
            for path in input_paths
        ]
        return [future.result() for future in futures]


def run_gulp(
    ncpus: int,
    input_path: Path,
    directory: Path,
    executable: Path | None,
    launcher: str,
    timeout: float | None = None,
) -> Path:
    """Peform a single GULP run."""
    launchers = {
        "mpi": ["mpirun", "-n", str(ncpus)],
        "srun": [
            "srun",
            "--exclusive",
            "--exact",
            "--nodes=1",
            "--ntasks",
            str(ncpus),
            "--cpus-per-task=1",
            "--distribution=block:block",
            "--hint=nomultithread",
        ],
    }
    l = "mpi" if launcher is None else launcher
    if l not in launchers:
        raise ValueError("Unable to find the specified launcher.")
    command = [*launchers[l], str(executable), input_path.stem]
    res = subprocess.run(
        command, cwd=directory, timeout=timeout, capture_output=True, check=True
    )
    gout_path = directory / f"{input_path.stem}.gout"

    return gout_path


def initialise_reference(
    settings: Settings,
    crystal: Structure,
    library_text: str,
) -> Structure:
    """Initialise the two reference calculation, return a Structure object containing relaxed Mott-Littleton Region I atoms."""

    master_input_path = settings.directory / "master.gin"
    reference_defect_path = settings.directory / "defect.gin"
    write_master_input(
        path=master_input_path,
        crystal=crystal,
        options=settings.options,
        library_text=library_text,
    )
    master_path = run_gulp(
        ncpus=settings.ncpus,
        input_path=master_input_path,
        directory=settings.directory,
        executable=settings.gulp_path,
        timeout=1000,
        launcher=settings.launcher,
    )

    output = read_output(master_path)
    final_crystal = Structure(cell=output.cell, atoms=output.atoms)

    write_ml_reference_input(
        path=reference_defect_path,
        crystal=final_crystal,
        settings=settings,
        library_text=library_text,
    )
    defect_dummy_out_path = run_gulp(
        ncpus=settings.ncpus,
        input_path=reference_defect_path,
        directory=settings.directory,
        executable=settings.gulp_path,
        timeout=1000,
        launcher=settings.launcher,
    )
    restart_path = defect_dummy_out_path.parent / "1.res"
    region1_atoms, centre_cart = read_restart(restart_path)

    reference_structure = Structure(
        cell=final_crystal.cell,
        atoms=final_crystal.atoms,
        region1_atoms=region1_atoms,
        centre=centre_cart,
    )
    return reference_structure
