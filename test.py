from interface.qe import read_qe
from structure.atoms import Atoms
from structure.cells import get_supercell, Primitive, print_cell
import numpy as np

symprec=1e-5
supercell_matrix = [[2, 0, 0],[0, 2, 0],[0, 0, 2]] # q vectors (integers)
primitive_matrix = [[1.0, 0.0, 0.0],[0.0, 1.0, 0.0], [0.0, 0.0, 1.0]] # primitive axis

unitcell = read_qe("scf.out") # returns Atoms object
#print unitcell
supercell = get_supercell(unitcell, supercell_matrix, symprec)
inv_supercell_matrix = np.linalg.inv(supercell_matrix)
primitive = Primitive(supercell, np.dot(inv_supercell_matrix, primitive_matrix), symprec)

self._modulation = Modulation(self._dynamical_matrix,
                                      self._primitive,
                                      dimension=dimension,
                                      phonon_modes=phonon_modes,
                                      factor=self._factor)


