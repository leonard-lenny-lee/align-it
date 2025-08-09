__all__ = [
    "prepare_mol", "prepare_mol_from_smiles", "decompose_mol", "recompose_mol"
]

from multiprocessing.pool import ThreadPool

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdRGroupDecomposition as rdRGD


def batch_prepare_mols_from_smiles(smiles: list[str], prepare_db: bool = False) -> list[Chem.Mol]:
    args = [(smi, prepare_db) for smi in smiles]
    with ThreadPool() as pool:
        mols = pool.starmap(prepare_mol_from_smiles, args)
    return mols


def prepare_mol_from_smiles(smi: str, prepare_db: bool = False) -> Chem.Mol:
    mol = Chem.MolFromSmiles(smi, sanitize=True)
    if mol is None:
        raise ValueError(f"Invalid SMILES {smi}.")
    mol = Chem.AddHs(mol)
    mol = prepare_mol(mol, prepare_db)
    return mol


def prepare_mol(mol: Chem.Mol, prepare_db: bool = False) -> Chem.Mol:
    # Replace wildcard atoms with protons for 3D embedding
    mol = Chem.RWMol(mol)
    wildcard_idxs = []

    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0 or (atom.GetAtomicNum() == 1 and prepare_db):
            atom.SetAtomicNum(1)
            wildcard_idxs.append(atom.GetIdx())

    # Embedding
    Chem.SanitizeMol(mol)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)

    # Replace protons back with wildcard atoms
    mol = Chem.RWMol(mol)

    for idx in wildcard_idxs:
        atom = mol.GetAtomWithIdx(idx)
        atom.SetAtomicNum(0)

    return mol.GetMol()


def decompose_mol(mol_smi: str, core_smi: str) -> tuple[Chem.Mol, list[Chem.Mol]]:
    mol = Chem.MolFromSmiles(mol_smi)
    core = Chem.MolFromSmiles(core_smi)

    if not mol.HasSubstructMatch(core):
        raise ValueError(f"{core_smi} is not a substructure of {mol_smi}")

    decomp = rdRGD.RGroupDecompose([core], [mol], asRows=False)[0]
    parts = [mol[0] for mol in decomp.values()]
    core, *r_groups = parts
    return core, r_groups


def recompose_mol(core: Chem.Mol, rgroups: list[Chem.Mol]) -> Chem.Mol:
    mol = core
    for a in mol.GetAtoms():
        if a.GetAtomicNum() == 0 and a.GetAtomMapNum() == 0:
            a.SetAtomicNum(1)
    if not rgroups:
        return mol
    combined = rgroups[0]
    for rgroup in rgroups[1:]:
        combined = AllChem.CombineMols(combined, rgroup)
    mol = AllChem.molzip(mol, combined)
    mol = AllChem.RemoveHs(mol)
    AllChem.Compute2DCoords(mol)
    return mol
