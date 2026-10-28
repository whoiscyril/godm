from collections import Counter
import numpy as np
import copy
import os

from ..helpers import check_collision, assign_charges, generate_shells
from pathlib import Path

from ..models import Settings, Species, Structure, Atom
from ..gulp.inputs import write_rq_input


def generate_structure(
    settings: Settings, reference: Structure, species: Species, seed: int, max_try: int
) -> Structure:
    """Generate a single structure after applying random quench procedures in Mott-Littleton Region I"""

    defect_list = settings.defects
    centre_defect = settings.centre_defect
    nimpure = Counter((d.label, d.host) for d in defect_list if d.kind == "impurity")
    ninter = Counter(d.label for d in defect_list if d.kind == "interstitial")
    nvac = Counter(d.label for d in defect_list if d.kind == "vacancy")
    r0 = settings.r0cutoff
    centre = np.array(reference.centre, dtype=float)

    rng = np.random.RandomState(seed)
    cores = [
        Atom(
            label=atom.label, kind=atom.kind, position=atom.position, charge=atom.charge
        )
        for atom in reference.region1_atoms
        if atom.kind == "core"
    ]

    positions = np.array([atom.position for atom in cores], dtype=float)
    distances = np.linalg.norm(positions - centre, axis=1)
    # Handle Impurity
    if centre_defect.kind == "impurity":
        key = (centre_defect.label, centre_defect.host)

        if key in nimpure and nimpure[key] > 0:
            nimpure[key] -= 1
    for (i, j), n in nimpure.items():
        if n == 0:
            continue
        eligible = [
            k for k, atom in enumerate(cores) if atom.label == j and distances[k] <= r0
        ]
        if len(eligible) < n:
            raise ValueError(
                f"Requested {n} {i} for {j} substitutions but only len({eligible}) were found"
            )

        for _ in range(n):
            index = int(rng.choice(eligible))
            cores[index].label = i
            eligible.remove(index)

    # Handle vacancies
    to_remove = set()
    if centre_defect.kind == "vacancy":
        key = centre_defect.label
        if key in nvac and nvac[key] > 0:
            nvac[key] -= 1
    for i, n in nvac.items():
        if n == 0:
            continue
        eligible = [
            k for k, atom in enumerate(cores) if atom.label == i and distances[k] <= r0
        ]

        if len(eligible) < n:
            raise ValueError(
                f"Requested {n} {i} vacancies but only len({eligible}) were found)"
            )
        for _ in range(n):
            index = int(rng.choice(eligible))
            to_remove.add(index)
            eligible.remove(index)
    cores = [atom for k, atom in enumerate(cores) if k not in to_remove]

    # Handle interstitial
    to_add = []
    if centre_defect.kind == "interstitial":
        key = centre_defect.label
        if key in ninter and ninter[key] > 0:
            ninter[key] -= 1
    for i, n in ninter.items():
        if n == 0:
            continue
        for _ in range(n):
            for attempt in range(max_try):
                r = r0 * rng.uniform(0, 1) ** (1 / 3)
                theta = np.arccos(rng.uniform(-1, 1))
                phi = rng.uniform(0, 2 * np.pi)
                d = np.array(
                    [
                        r * np.sin(theta) * np.cos(phi),
                        r * np.sin(theta) * np.sin(phi),
                        r * np.cos(theta),
                    ]
                )
                pos = centre + d
                if np.linalg.norm(centre - pos) <= r0:
                    if any(
                        check_collision(pos, atom.position, 1.0)
                        for atom in cores + to_add
                    ):
                        continue

                    x, y, z = map(float, pos)
                    a = Atom(label=i, kind="core", position=(x, y, z))
                    to_add.append(a)
                    break
            else:
                raise RuntimeError(
                    f"Unable to place interstitial {i} after {max_try} attempts"
                )
    cores.extend(to_add)

    # Add shells
    cs = generate_shells(cores, species)

    # Assign Charges
    cs = assign_charges(cs, species)

    shell_pairs = [
        (index - 1, index) for index, atom in enumerate(cs) if atom.kind == "shell"
    ]

    return Structure(
        cell=reference.cell,
        atoms=copy.deepcopy(reference.atoms),
        region1_atoms=cs,
        centre=reference.centre,
        shell_pairs=shell_pairs,
    )


def generate_rq_batch(
    settings: Settings,
    reference: Structure,
    library_text: str,
    species: Species,
) -> list[Path]:
    """Generate and write multiple Mott-Littleton random quench structures"""
    batch_rng = np.random.RandomState(17)
    seeds = batch_rng.randint(0, 2**32, size=settings.rq_steps, dtype=np.uint32)
    run_dir = settings.directory / "run"
    sim_file_paths = []
    try:
        os.makedirs(run_dir)
    except FileExistsError as e:
        raise FileExistsError(
            "run directory exists, consider renaming it or delet it"
        ) from e
    for index, seed in enumerate(seeds):
        candidate = generate_structure(
            settings=settings,
            reference=reference,
            species=species,
            seed=int(seed),
            max_try=1000,
        )
        filename = f"{index}.gin"
        file_dir = run_dir / filename.removesuffix(".gin")
        file_dir.mkdir(parents=True, exist_ok=False)
        file_path = file_dir / filename
        sim_file_paths.append(file_path)
        write_rq_input(
            path=file_path,
            structure=candidate,
            settings=settings,
            library_text=library_text,
            seed=seed,
        )
    return sim_file_paths
