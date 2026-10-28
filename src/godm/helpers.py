import numpy as np

from .models import Atom, Species

CORE_REP = ["cor", "c", "core"]
SHELL_REP = ["shel", "she", "s", "shell"]


def check_collision(
    pos1: tuple[float, float, float], pos2: tuple[float, float, float], cutoff: float
) -> bool:
    """Check if positions of two atoms are too close given a cutoff"""
    a1 = np.array(pos1)
    a2 = np.array(pos2)
    return np.linalg.norm(a1 - a2) <= cutoff


def assign_charges(atoms: list[Atom], species: Species) -> list[Atom]:
    """Assign charges to atoms/ions cores and shells according to the library file"""
    for atom in atoms:
        key = (atom.label, atom.kind)
        if key not in species:
            raise ValueError(f"Charge for species {key} not found")
        atom.charge = species[key]
    return atoms


def generate_shells(atoms: list[Atom], species: Species) -> list[Atom]:
    """Generate corresponding shells for atoms/ions if they use the shell model"""
    result = []

    for atom in atoms:
        if atom.kind != "core":
            raise ValueError("generate_shells expects core atoms only")

        result.append(atom)
        key = (atom.label, "shell")
        if key in species:
            shell = Atom(
                label=atom.label,
                kind="shell",
                position=atom.position,
                charge=species[key],
            )
            result.append(shell)

    return result


def normalise_core_shell(kind: str):
    """Normalise different core/shell representations in GULP inputs and outputs"""
    kind = kind.lower()
    if kind in CORE_REP:
        return "core"
    if kind in SHELL_REP:
        return "shell"
    raise ValueError("Unknown atom type: check core/shell representation")
