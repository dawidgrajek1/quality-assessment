from Bio import PDB
from Bio.PDB.PDBIO import PDBIO

parser = PDB.PDBParser(QUIET=True)

model = parser.get_structure("model", "zad2output.pdb")
ref = parser.get_structure("ref", "R1107_reference.pdb")

for chain1, chain2 in zip(model.get_chains(), ref.get_chains()):
    chain1.id = chain2.id

for a1, a2 in zip(model.get_atoms(), ref.get_atoms()):
    a1.set_serial_number(a2.get_serial_number())

io = PDBIO()
io.set_structure(model)
io.save("modelFinal.pdb")
io.set_structure(ref)
io.save("refFinal.pdb")
