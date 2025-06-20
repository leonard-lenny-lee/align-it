__all__ = ["view_pharmacophores"]

import pathlib
import subprocess

import numpy as np
from rdkit import Chem
import pymol2
from pymol import cgo

from .cpyalignit import *
from .draw import PHARM_COLORS

DEFAULT_SPHERE_ALPHA = 0.5
DEFAULT_SPHERE_RADIUS = 1
DEFAULT_ARROW_RADIUS = 0.1

type Model = tuple[str, Pharmacophore | None, Chem.Mol | None]
type Vec3 = tuple[float, float, float]


def view_pharmacophores(
    session_name: str,
    models: list[Model],
    launch_pymol: bool = True,
) -> None:
    """View models (A pharmacophore and mol) in Pymol"""
    path = pathlib.Path(session_name)
    pse_path = str(path.with_suffix(".pse"))

    with pymol2.PyMOL() as pymol:
        for model in models:
            _view_pharmacophore(session_name, model, pymol)

        pymol.cmd.set("transparency_mode", 1)
        pymol.cmd.save(pse_path)

    if launch_pymol:
        subprocess.run(["pymol", pse_path])


def _view_pharmacophore(
    session_name: str, model: Model, instance: pymol2.PyMOL
) -> None:
    model_name, pharm, mol = model
    group_name = model_name
    instance.cmd.group(group_name)

    if mol is not None:
        pdb_name = f"{session_name}_{model_name}"
        pdb_path = pathlib.Path(pdb_name).with_suffix(".pdb")
        Chem.MolToPDBFile(mol, pdb_path)
        instance.cmd.load(pdb_path, pdb_name)
        instance.cmd.group(group_name, pdb_name, "add")

    if pharm is None:
        return

    for i, p in enumerate(pharm):
        cgo_list = []
        pos = p.point.x, p.point.y, p.point.z
        color = PHARM_COLORS[p.func]
        func_label = f"{model_name}_{p.func}_{i}"
        sphere = _sphere(pos, DEFAULT_SPHERE_RADIUS, color)
        cgo_list += sphere

        if p.hasNormal:
            normal_pos = p.normal.x, p.normal.y, p.normal.z
            reverse = p.func == FuncGroup.HACC
            arrow = _arrow(pos, normal_pos, DEFAULT_ARROW_RADIUS, color, reverse)
            cgo_list += arrow

        instance.cmd.load_cgo(cgo_list, func_label)
        instance.cmd.group(group_name, func_label, "add")


def _sphere(pos: Vec3, radius: float, hex: str) -> list[float]:
    rgb = _hex_to_rgb(hex)
    sphere = [
        cgo.ALPHA,
        DEFAULT_SPHERE_ALPHA,
        cgo.COLOR,
        *rgb,
        cgo.SPHERE,
        *pos,
        radius,
    ]
    return sphere


def _arrow(
    pos0: Vec3, pos1: Vec3, radius: float, hex: str, reverse: bool
) -> list[float]:
    rgb = _hex_to_rgb(hex)
    _pos0, _pos1 = np.array(pos0), np.array(pos1)
    _dir = _pos1 - _pos0
    _pos2 = _pos1 + _dir
    if reverse:
        _pos0, _pos2 = _pos2, pos0
    cylinder = [cgo.CYLINDER, *_pos0, *_pos1, radius, *rgb, *rgb]
    cone = [cgo.CONE, *_pos1, *_pos2, DEFAULT_ARROW_RADIUS * 2, 0, *rgb, *rgb, 1, 1]
    arrow = cylinder + cone
    return arrow


def _hex_to_rgb(hex: str) -> tuple[float, ...]:
    hex = hex.lstrip("#")
    return tuple(int(hex[i : i + 2], 16) / 255.0 for i in (0, 2, 4))
