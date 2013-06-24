# striped interface for QE
import re
import numpy as np
from modulation.units import Bohr
from structure.atoms import Atoms

def read_qe(filename, symbols=None):
    alat = 0.0
    cell = []
    output = open(filename)
    while True:
        line = output.readline()
        if not line:
            break

        p = re.search(r'lattice parameter \(alat\)\s+=\s+([\d\.]+)', line)
        if p:
            alat = float(p.group(1))
            print "Found alat %f" % alat
            continue # while True

        p = re.search(r'crystal axes: \(cart\. coord\. in units of alat\)', line)
        if p:
            for i in range(0, 3):
                p1 = re.search(r'.+?(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)', output.readline())
                cell.append([ float(p1.group(1)), float(p1.group(2)), float(p1.group(3)) ])
            cell = np.array(cell)# * scale
            print "Crystal axes in units of alat:"
            print cell

            cell = cell*alat*Bohr # [A] now
            continue # while True

        p = re.search(r'Crystallographic axes', line)
        if p:
            output.readline()
            output.readline()

            positions = []
            expanded_symbols = []
            while True:
                p1 = re.search(r'^\s*\d+\s+(\w+).+?(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)', output.readline())
                if not p1:
                    break
                expanded_symbols.append( p1.group(1).upper() )
                positions.append([ float(p1.group(2)), float(p1.group(3)), float(p1.group(4)) ])

            atoms = Atoms(symbols=expanded_symbols, cell=cell, positions=positions)
            return atoms
