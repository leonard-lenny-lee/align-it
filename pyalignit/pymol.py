import pathlib
import subprocess

import numpy as np
from rdkit import Chem
import pymol2
from pymol import cgo

from .cpyalignit import *
from .draw import PHARM_COLORS

SPHERE_ALPHA = 0.5
SPHERE_RADIUS = 1
ARROW_RADIUS = 0.1

type Model = tuple[Pharmacophore, Chem.Mol | None]
type _Coord = tuple[float, float, float]


def view_pharmacophore(
    name: str,
    pharm: Pharmacophore,
    mol: Chem.Mol | None = None,
    launch_pymol: bool = False,
) -> None:
    path = pathlib.Path(name)
    pdb_path = str(path.with_suffix(".pdb"))
    pse_path = str(path.with_suffix(".pse"))

    if mol is not None:
        Chem.MolToPDBFile(mol, pdb_path)

    with pymol2.PyMOL() as pymol:
        pymol.cmd.load(pdb_path)

        for i, p in enumerate(pharm):
            cgo_list = []
            pos = p.point.x, p.point.y, p.point.z
            color = PHARM_COLORS[p.func]
            label = f"{i+1}-{p.func}"
            sphere = _sphere(pos, SPHERE_RADIUS, color)
            cgo_list += sphere

            if p.hasNormal:
                normal_pos = p.normal.x, p.normal.y, p.normal.z
                reverse = p.func == FuncGroup.HACC
                arrow = _arrow(pos, normal_pos, ARROW_RADIUS, color, reverse)
                cgo_list += arrow

            pymol.cmd.load_cgo(cgo_list, f"{label}")

        pymol.cmd.set("transparency_mode", 1)
        pymol.cmd.save(pse_path)

    if launch_pymol:
        subprocess.run(["pymol", pse_path])


def _sphere(pos: _Coord, radius: float, hex: str) -> list[float]:
    rgb = _hex_to_rgb(hex)
    sphere = [cgo.ALPHA, SPHERE_ALPHA, cgo.COLOR, *rgb, cgo.SPHERE, *pos, radius]
    return sphere


def _arrow(
    pos0: _Coord, pos1: _Coord, radius: float, hex: str, reverse: bool
) -> list[float]:
    rgb = _hex_to_rgb(hex)
    _pos0, _pos1 = np.array(pos0), np.array(pos1)
    _dir = _pos1 - _pos0
    _pos2 = _pos1 + _dir
    if reverse:
        _pos0, _pos2 = _pos2, pos0
    cylinder = [cgo.CYLINDER, *_pos0, *_pos1, radius, *rgb, *rgb]
    cone = [cgo.CONE, *_pos1, *_pos2, ARROW_RADIUS * 2, 0, *rgb, *rgb, 1, 1]
    arrow = cylinder + cone
    return arrow


def _hex_to_rgb(hex: str) -> tuple[float, ...]:
    hex = hex.lstrip("#")
    return tuple(int(hex[i : i + 2], 16) / 255.0 for i in (0, 2, 4))
