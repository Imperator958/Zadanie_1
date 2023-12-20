import Bio.PDB
import numpy as np
from Bio.SVDSuperimposer import SVDSuperimposer

def rmsd(coords_1, coords_2):
    diffrence = coords_1 - coords_2
    return np.sqrt(np.mean(np.sum(diffrence**2, axis=1)))

pdb_parser = Bio.PDB.PDBParser(QUIET=True)

sample = pdb_parser.get_structure("sample", "R1107TS081.pdb")
ref = pdb_parser.get_structure("reference", "R1107_reference.pdb")

ref_atoms = []
sample_atoms = []

get_ref_atoms = ref.get_atoms()
get_sample_atoms = sample.get_atoms()

for a in get_ref_atoms:
    ref_atoms.append(a)
for a in get_sample_atoms:
    sample_atoms.append(a)

ref_models = [ref_atoms[0:1461], ref_atoms[0:1464], ref_atoms[0:1461], ref_atoms[0:1461]]
sample_models = [sample_atoms[3:1464], sample_atoms[1467:2928], sample_atoms[2931:4392], sample_atoms[4395:5856]]

ref_model = ref[0]
sample_model = sample[0]

i = 0
while i < 4:

    super_imposer = Bio.PDB.Superimposer()
    super_imposer.set_atoms(ref_models[i], sample_models[i])
    super_imposer.apply(sample_model.get_atoms())

    coords_ref = np.array(([atom.get_coord() for atom in ref_models[i]]))
    coords_sample = np.array(([atom.get_coord() for atom in sample_models[i]]))

    sup = SVDSuperimposer()
    sup.set(coords_ref, coords_sample)
    sup.run()
    ref_coords = sup.reference_coords
    transformed = sup.get_transformed()
    get_rmsd = rmsd(transformed, ref_coords)
    print("RMSD dla modelu {0}: ".format(i+1), get_rmsd)

    i += 1

