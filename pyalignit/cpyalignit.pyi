from enum import StrEnum
from typing import Any

def GetVersion() -> str: ...

class Coordinate:
    x: float
    y: float
    z: float

class FuncGroup(StrEnum):
    AROM = "AROM"
    HDON = "HDON"
    HACC = "HACC"
    LIPO = "LIPO"
    POSC = "POSC"
    NEGC = "NEGC"
    HYBH = "HYBH"
    HYBL = "HYBL"
    EXCL = "EXCL"
    UNDEF = "UNDEF"
    ATTA = "ATTA"

class PharmacophorePoint:
    point: Coordinate
    normal: Coordinate
    func: FuncGroup
    alpha: float
    hasNormal: bool

type Pharmacophore = list[PharmacophorePoint]

def CalcPharmacophore(
    mol: Any,
    calcArom: bool = True,
    calcHDon: bool = True,
    calcHAcc: bool = True,
    calcLipo: bool = True,
    calcCharge: bool = True,
    calcHybrid: bool = True,
) -> Pharmacophore: ...
def MergePharmacophore(pharmacophore: Pharmacophore) -> Pharmacophore: ...
def AlignPharmacophore(
    ref: Pharmacophore,
    probe: Pharmacophore,
    probeMol: Any = None,
    epsilon: float = 0.5,
    useNormals: bool = True,
    useExclusion: bool = False,
) -> tuple[float, float, float]: ...
def AlignMol(
    ref: Any,
    probe: Any,
    calcArom: bool = True,
    calcHDon: bool = True,
    calcHAcc: bool = True,
    calcLipo: bool = True,
    calcCharge: bool = True,
    calcHybrid: bool = True,
    merge: bool = False,
    epsilon: float = 0.5,
    useNormals: bool = True,
    useExclusion: bool = False,
) -> tuple[float, float, float]: ...
