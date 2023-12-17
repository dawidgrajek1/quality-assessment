from Bio import PDB
from sys import argv

atom_radii = {
    "H": 1.20,
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "S": 1.80,
    "F": 1.47,
    "P": 1.80,
    "CL": 1.75,
    "MG": 1.73,
}


def count_clashes(structure, clash_cutoff=0.4):
    counter = 0
    clashCutoffs = {
        i + "_" + j: (atom_radii[i] + atom_radii[j] - clash_cutoff)
        for i in atom_radii
        for j in atom_radii
    }

    atoms = [x for x in structure.get_atoms() if x.element in atom_radii]
    print("Number of atoms: ", len(atoms))
    for i in range(0, len(atoms)):
        for j in range(i + 1, len(atoms)):
            if atoms[i].parent.id == atoms[j].parent.id:
                continue

            # elif (atoms[j].element == "C" and atoms[i].element == "P") or (
            #     atoms[j].element == "P" and atoms[i].element == "C"
            # ):
            #     continue

            # elif (atoms[j].element == "C" and atoms[i].element == "P") or (
            #     atoms[j].element == "P" and atoms[i].element == "C"
            # ):
            #     continue

            # elif (atoms[j].element == "C" and atoms[i].element == "O") or (
            #     atoms[j].element == "O" and atoms[i].element == "C"
            # ):
            #     continue

            if (
                atoms[i] - atoms[j]
                <= clashCutoffs[atoms[i].element + "_" + atoms[j].element]
            ):
                # print(
                #     atoms[i].serial_number,
                #     atoms[i].element,
                #     "\t",
                #     atoms[j].serial_number,
                #     atoms[j].element,
                #     "\t",
                #     clashCutoffs[atoms[i].element + "_" + atoms[j].element]
                #     - (atoms[i] - atoms[j]),
                # )
                counter += 1
    print("Number of clashes: ", counter)
    return counter * 1000 / len(atoms)


if len(argv) != 2:
    print("Usage: python clashScore.py <pdb file>")
    exit()

parser = PDB.PDBParser(QUIET=True)
structure = parser.get_structure("struct", argv[1])
print("Clash score: ", count_clashes(structure))
