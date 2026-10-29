import numpy as np
from slime_proteo_md.analysis import calculate_radius_of_gyration, kabsch_rmsd


def test_radius_of_gyration():
    coords = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]], dtype=float)
    rg = calculate_radius_of_gyration(coords)
    assert rg > 0
    assert rg < 1


def test_kabsch_rmsd_identical():
    coords = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]], dtype=float)
    assert kabsch_rmsd(coords, coords) < 1e-10


def test_kabsch_rmsd_rotated():
    p = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]], dtype=float)
    q = np.array([[0, 0, 0], [0, 1, 0], [-1, 0, 0]], dtype=float)
    assert kabsch_rmsd(p, q) < 1e-10
