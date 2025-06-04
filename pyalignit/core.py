from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdRGroupDecomposition as rdRGD


def prepare_mol(mol: Chem.Mol) -> Chem.Mol:
    # Replace dummy atoms with protons
    mol = Chem.RWMol(mol)
    dummy_atom_idxs = []

    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0:
            atom.SetAtomicNum(1)
            dummy_atom_idxs.append(atom.GetIdx())

    # Energy minimization
    Chem.SanitizeMol(mol)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    AllChem.MMFFOptimizeMolecule(mol)

    # Replace protons back with dummy atoms
    mol = Chem.RWMol(mol)

    for idx in dummy_atom_idxs:
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
