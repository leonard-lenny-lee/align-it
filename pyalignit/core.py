from rdkit import Chem
from rdkit.Chem import AllChem


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
