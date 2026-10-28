import shlex
from pathlib import Path

from .models import Atom, Defect, Settings, Structure


def read_inputs() -> tuple[Settings, Structure]:
    "Read input/CONTROL and input/CONFIG and return settings plus the crystal."

    p = Path(".")
    q = p / "inputs"
    required = ["CONTROL", "CONFIG"]
    _control_file = q / required[0]
    _config_file = q / required[1]

    control_vals = _parse_control(_control_file)
    config_vals = _parse_config(_config_file)
    settings_vals = {
        **control_vals,
        "directory": p.resolve(),
        "options": config_vals["options"],
        "keywords": config_vals["keywords"],
    }

    settings = Settings(**settings_vals)
    structure = Structure(cell=config_vals["cell"], atoms=config_vals["atoms"])
    return (settings, structure)


def _parse_config(config_path: Path) -> dict:
    "Parse the CONFIG file"
    vals = {}
    in_kw = False
    in_cell = False
    in_atoms = False
    in_options = False
    options = []
    atoms = []
    with config_path.open(mode="r") as f:
        lines = f.readlines()
        for line in lines:
            tokens = line.split()
            if not tokens:
                continue
            if tokens[0] == "BEGIN_KW":
                in_kw = True
                continue
            if tokens[0] == "BEGIN_CELL":
                in_cell = True
                continue
            if tokens[0] == "BEGIN_ATOMS_FRAC":
                in_atoms = True
                continue
            if tokens[0] == "BEGIN_OPTIONS":
                in_options = True
                continue

            elif tokens[0] == "END_KW":
                in_kw = False
                continue
            elif tokens[0] == "END_CELL":
                in_cell = False
                continue

            elif tokens[0] == "END_ATOMS_FRAC":
                in_atoms = False
                continue
            elif tokens[0] == "END_OPTIONS":
                in_options = False
                continue

            elif in_kw:
                vals["keywords"] = tokens
            elif in_cell:
                a, b, c, alpha, beta, gamma = map(float, tokens)
                vals["cell"] = (a, b, c, alpha, beta, gamma)

            elif in_options:
                options.append(line)
            elif in_atoms:
                x, y, z = map(float, tokens[2:5])
                at = Atom(
                    label=tokens[0],
                    kind=normalise_core_shell(tokens[1]),
                    position=(x, y, z),
                )
                atoms.append(at)

        vals["options"] = options
        vals["atoms"] = atoms

    return vals


def _parse_control(
    control_path: Path,
) -> dict:
    vals = {}
    "Parse the CONTROL file"
    defects: list[Defect] = []
    in_dlist = False
    input_dir = control_path.resolve().parent
    with control_path.open(mode="r") as f:
        lines = f.readlines()
    for line in lines:
        tokens = line.split()

        if not tokens:
            continue

        if tokens[0] == "BEGIN_DLIST":
            in_dlist = True
            continue
        elif tokens[0] == "END_DLIST":
            in_dlist = False
            continue
        elif in_dlist:
            kind = tokens[0]
            if kind == "impurity" and len(tokens) == 3:
                defects.append(Defect(kind, label=tokens[1], host=tokens[2]))
            elif kind in {"interstitial", "vacancy"} and len(tokens) == 2:
                defects.append(Defect(kind, label=tokens[1]))
        if tokens[0] == "LAUNCHER":
            vals["launcher"] = tokens[-1].lower()
        if tokens[0] == "NUM_WORKERS":
            vals["num_workers"] = int(tokens[-1])
        if tokens[0] == "NCPUS":
            vals["ncpus"] = int(tokens[-1])

        if tokens[0] == "ML_SIZE":
            try:
                r1 = float(tokens[-2])
                r2 = float(tokens[-1])
            except ValueError:
                raise ValueError("Error parsing Mott-Littleton Region cut-offs.")

            vals["mlcutoffs"] = (r1, r2)

        elif tokens[0] == "R0_CUTOFF":
            try:
                cutoff = float(tokens[-1])
            except ValueError:
                raise ValueError("Error parsing Region 0 cut off value.")

            vals["r0cutoff"] = cutoff

        elif tokens[0] == "ALGORITHM":
            algo = tokens[-1]
            vals["algorithm"] = algo
        elif tokens[0] == "RQ_STEPS":
            try:
                steps = int(tokens[-1])
            except ValueError:
                raise ValueError("Error parsing random quench steps")
            vals["rq_steps"] = steps
        elif tokens[0] == "GULP_EXE":
            _, p = shlex.split(line)
            vals["gulp_path"] = Path(p).resolve()
        elif tokens[0] == "LIBRARY":
            _, p = shlex.split(line)
            vals["library_path"] = (input_dir / p).resolve()
        elif tokens[0] == "CENTRE_FRACTIONAL":
            try:
                x, y, z = float(tokens[-3]), float(tokens[-2]), float(tokens[-1])
            except ValueError:
                raise ValueError("Error parsing the defect centre coordinates")
            vals["centre_fractional"] = (x, y, z)

        elif tokens[0] == "CENTRE_DEFECT":
            if len(tokens) == 4:
                d = Defect(kind=tokens[1], label=tokens[2], host=tokens[3])
                vals["centre_defect"] = d
            elif len(tokens) == 3:
                d = Defect(kind=tokens[1], label=tokens[2])
                vals["centre_defect"] = d

    vals["defects"] = defects

    return vals
