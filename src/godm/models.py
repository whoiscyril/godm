from pathlib import Path

Species = dict[tuple[str, str], float]
Position = tuple[float, float, float]
Cell = tuple[float, float, float, float, float, float]


class Atom:
    def __init__(
        self,
        label: str,
        kind: str,
        position: Position,
        charge: float | None = None,
    ) -> None:
        self.label = label
        self.kind = kind
        self.position = position
        self.charge = charge

    def __str__(self):
        return f"{self.label}, {self.kind}, {self.position}, {self.charge}"


class Defect:
    def __init__(self, kind: str, label: str, host: str | None = None) -> None:
        self.kind = kind
        self.label = label
        self.host = host


class Settings:
    def __init__(
        self,
        *,
        directory: Path,
        algorithm: str,
        rq_steps: int,
        gulp_path: Path | None,
        library_path: Path,
        num_workers: int,
        ncpus: int,
        launcher: str,
        mlcutoffs: tuple[float, float],
        r0cutoff: float,
        centre_fractional: Position,
        centre_defect: Defect,
        defects: list[Defect],
        options: list[str],
        keywords: list[str],
    ) -> None:

        self.directory = directory
        self.algorithm = algorithm
        self.rq_steps = rq_steps
        self.gulp_path = gulp_path
        self.library_path = library_path
        self.mlcutoffs = mlcutoffs
        self.r0cutoff = r0cutoff
        self.centre_fractional = centre_fractional
        self.centre_defect = centre_defect
        self.defects = defects
        self.options = options
        self.keywords = keywords
        self.launcher = launcher
        self.ncpus = ncpus
        self.num_workers = num_workers

    def __str__(self) -> str:
        options = "\n".join("    " + line.rstrip("\r\n") for line in self.options)
        return (
            "Settings\n"
            f"  Directory: {self.directory}\n"
            f"  Algorithm: {self.algorithm}\n"
            f"  RQ steps: {self.rq_steps}\n"
            f"  GULP executable: {self.gulp_path}\n"
            f"  Library: {self.library_path}\n"
            f"  ML cutoffs: {self.mlcutoffs}\n"
            f"  Region 0 cutoff: {self.r0cutoff}\n"
            f"  Fractional centre: {self.centre_fractional}\n"
            f"  Centre defect: {self.centre_defect.kind} {self.centre_defect.label}"
            f"{(' ' + self.centre_defect.host) if self.centre_defect.host else ''}\n"
            f"  Defects: {len(self.defects)}\n"
            f"  Keywords: {' '.join(self.keywords)}\n"
            f"  Options:\n{options}"
        )


class Structure:
    def __init__(
        self,
        cell: Cell,
        atoms: list[Atom],
        region1_atoms: list[Atom] | None = None,
        centre: Position | None = None,
        shell_pairs: list[tuple[int, int]] | None = None,
    ) -> None:

        self.cell = cell
        self.atoms = atoms
        self.region1_atoms = region1_atoms
        self.centre = centre
        self.shell_pairs = shell_pairs

    def __str__(self) -> str:
        atom_list = "\n".join(f"{atom}" for atom in self.atoms)
        return (
            "Structure\n"
            f"  Cell (a, b, c, alpha, beta, gamma): {self.cell}\n"
            f"  Atom coordinates: \n"
            f"  {atom_list}\n"
            f"  Region 1 atoms: {len(self.region1_atoms) if self.region1_atoms is not None else 0}\n"
            f"  Cartesian centre: {self.centre}\n"
            f"  Shell pairs: {len(self.shell_pairs) if self.shell_pairs is not None else 0}"
        )


class GulpOutput:
    def __init__(
        self,
        *,
        finished: bool = False,
        converged: bool | None = None,
        has_bad_energy: bool = False,
        lattice_energy: float | None = None,
        defect_energy: float | None = None,
        defect_gnorm: float | None = None,
        lattice_gnorm: float | None = None,
        cpu_time: float | None = None,
        cell: Cell | None = None,
        atoms: list[Atom] | None = None,
        region1_atoms: list[Atom] | None = None,
        centre: Position | None = None,
    ) -> None:

        self.finished = finished
        self.converged = converged
        self.has_bad_energy = has_bad_energy
        self.lattice_energy = lattice_energy
        self.defect_energy = defect_energy
        self.lattice_gnorm = lattice_gnorm
        self.defect_gnorm = defect_gnorm
        self.cpu_time = cpu_time
        self.cell = cell
        self.atoms = atoms
        self.region1_atoms = region1_atoms
        self.centre = centre
