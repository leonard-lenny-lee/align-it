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
    EXIT = "EXIT"

class PharmacophorePoint:
    point: Coordinate
    normal: Coordinate
    func: FuncGroup
    alpha: float
    hasNormal: bool

class Result:
    refId: str
    refVolume: float
    dbId: str
    dbVolume: float
    overlapVolume: float
    exclVolume: float
    resPharSize: int
    tanimoto: float
    tversky_ref: float
    tversky_db: float
    rankbyScore: float
    resMol: Any
    resPhar: Pharmacophore

Pharmacophore = list[PharmacophorePoint]

def CalcPharmacophore(
    mol: Any,
    calcArom: bool = True,
    calcHDon: bool = True,
    calcHAcc: bool = True,
    calcLipo: bool = True,
    calcCharge: bool = True,
    calcHybrid: bool = True,
    calcExits: bool = True,
) -> Pharmacophore: ...
def MergePharmacophore(pharmacophore: Pharmacophore) -> Pharmacophore: ...
def AlignPharmacophore(
    ref: Pharmacophore,
    probe: Pharmacophore,
    probeMol: Any = None,
    epsilon: float = 0.5,
    useNormals: bool = True,
    useExclusion: bool = False,
) -> Result: ...
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
    calcExits: bool = True,
) -> Result: ...
