from pathlib import Path

from ..models import Settings, Species, Structure
from ..helpers import normalise_core_shell


def read_library(path: Path) -> tuple[str, Species]:
    """Read the potential library and return the original library as a string and a species-charge dictionary."""
    with path.open(mode="r") as f:
        content = f.readlines()
        to_parse = 0
        s: Species = {}
    for line in content:
        tokens = line.split()
        if not tokens:
            continue
        if tokens[0] == "species":
            to_parse = int(tokens[-1])
            continue
        if to_parse > 0:
            label, kind, charge = tokens
            kind = normalise_core_shell(kind)
            s[(label, kind)] = float(charge)
            to_parse -= 1
    c = "".join(content)
    return (c, s)


def write_master_input(
    path: Path, crystal: Structure, options: list[str], library_text: str
) -> None:
    """Write the .gin file for master bulk relaxation"""
    with path.open(mode="w", encoding="utf-8") as f:
        f.write("opti conp nosymm \n")
        f.write("cell\n")
        f.write(" ".join(str(s) for s in crystal.cell) + "\n")
        f.write("fractional \n")
        for atom in crystal.atoms:
            c = " ".join(str(s) for s in atom.position)
            line = f"{atom.label} {atom.kind} {c}"
            f.write(line + "\n")
        for option in options:
            f.write(option.rstrip("\r\n") + "\n")
        f.write(library_text)


def write_ml_reference_input(
    path: Path, crystal: Structure, settings: Settings, library_text: str
) -> None:
    """Write the dummy defect calculation used to obtain the reference Mott-Littleton Region I."""
    with path.open(mode="w", encoding="utf-8") as f:
        f.write("defect nodsymm energy regi_before\n")
        f.write("cell\n")
        f.write(" ".join(str(s) for s in crystal.cell) + "\n")
        f.write("fractional \n")
        for atom in crystal.atoms:
            c = " ".join(str(s) for s in atom.position)
            line = f"{atom.label} {atom.kind} {c}"
            f.write(line + "\n")
        f.write(f"size {settings.mlcutoffs[0]} {settings.mlcutoffs[1]} \n")
        defect_coord = " ".join(str(s) for s in settings.centre_fractional)
        f.write(
            f"{settings.centre_defect.kind} {settings.centre_defect.label} {defect_coord} \n"
        )
        f.write(f"centre {defect_coord}\n")
        for option in settings.options:
            f.write(option.rstrip("\r\n") + "\n")
        f.write(library_text)


def write_rq_input(
    path: Path, structure: Structure, settings: Settings, library_text: str, seed: int
) -> None:
    """Write a single input file after random quench / basin-hopping moves have been proposed."""
    with path.open(mode="w", encoding="utf-8") as f:
        f.write(f"#seed {seed}\n")
        f.write("opti conp defect nodsymm freq\n")
        f.write("cell\n")
        f.write(" ".join(str(s) for s in structure.cell) + "\n")
        f.write("fractional \n")
        for atom in structure.atoms:
            c = " ".join(str(s) for s in atom.position)
            line = f"{atom.label} {atom.kind} {c}"
            f.write(line + "\n")
        f.write(f"centre cart {' '.join(str(s) for s in structure.centre)}\n")
        f.write(f"size {settings.mlcutoffs[0]} {settings.mlcutoffs[1]}\n")
        f.write("region_1\n")
        for atom in structure.region1_atoms:
            line = f"{atom.label} {atom.kind} {' '.join(str(s) for s in atom.position)}"
            f.write(line + "\n")
        for option in settings.options:
            f.write(option.rstrip("\r\n") + "\n")
        f.write(library_text)
