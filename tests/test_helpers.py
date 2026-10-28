from godm.helpers import check_collision, normalise_core_shell
import pytest


@pytest.mark.parametrize(
    "pos1, pos2, cutoff, expected",
    [
        ((0, 0, 0), (0, 0, 0), 5.0, True),
        ((0, 0, 0), (3, 0, 0), 5.0, True),
        ((0, 0, 0), (3, 4, 0), 5.0, True),
        ((0, 0, 0), (6, 0, 0), 5.0, False),
    ],
)
def test_check_collisions(pos1, pos2, cutoff, expected):
    assert check_collision(pos1, pos2, cutoff) == expected


@pytest.mark.parametrize("foo, expected", [("cor", "core"), ("s", "shell")])
def test_normalise_core_shell(foo, expected):
    assert normalise_core_shell(foo) == expected


def test_normalise_core_shell_exception():
    with pytest.raises(ValueError, match="Unknown atom type"):
        normalise_core_shell("cs")
