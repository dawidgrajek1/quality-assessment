from Bio import PDB
from Bio.PDB.Superimposer import Superimposer
from math import sqrt


def calculateRMSD(A: list, B: list) -> float:
    sum: float = 0

    for a, b in zip(A, B):
        sum += (a - b) ** 2
    return sqrt(sum / len(A))


parser = PDB.PDBParser(QUIET=True)
si = Superimposer()

ref = parser.get_structure("ref", "R1107_reference.pdb")
refAtoms = [x for x in ref.get_atoms()]
models = parser.get_structure("models", "R1107TS081.pdb")

print("Model\tRMSD")
for model in models:
    modelAtoms = [x for x in model.get_atoms()]
    modelAtoms = modelAtoms[abs(len(modelAtoms) - len(refAtoms)):]

    si.set_atoms(refAtoms, modelAtoms)
    si.apply(modelAtoms)

    print(f"{str(model)[-2]}\t{calculateRMSD(refAtoms, modelAtoms)}")
