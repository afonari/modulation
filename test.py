from interface.qe import read_qe, read_qe_dynmat
from structure.atoms import Atoms
from structure.cells import get_supercell, Primitive, print_cell
from modulation import Modulation
import numpy as np

symprec=1e-5
supercell_matrix = [[6, 0, 0],[0, 6, 0],[0, 0, 6]] # q vectors (integers)
primitive_matrix = [[1.0, 0.0, 0.0],[0.0, 1.0, 0.0], [0.0, 0.0, 1.0]] # primitive axis

unitcell, alat = read_qe("scf.out") # returns Atoms object
natoms = len(unitcell.get_scaled_positions())
supercell = get_supercell(unitcell, supercell_matrix, symprec)
inv_supercell_matrix = np.linalg.inv(supercell_matrix)
primitive = Primitive(supercell, np.dot(inv_supercell_matrix, primitive_matrix), symprec)

eigval, eigenvec, q = read_qe_dynmat("dynmat.dyn2", natoms, alat)
q = [0.0, 0.0, 0.1666667]
print eigval
print eigenvec
#mod = Modulation(q, eigval, eigenvec, primitive, [1, 1, 1])
#mod.write_yaml()
#mod.write()


