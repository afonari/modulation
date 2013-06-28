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
                p1 = re.search(r'^\s*(\w+)\s*(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)', output.readline())
                if not p1:
                    break
                expanded_symbols.append( p1.group(1) )
                positions.append([ float(p1.group(2)), float(p1.group(3)), float(p1.group(4)) ])

            atoms = Atoms(symbols=expanded_symbols, cell=cell, positions=positions)
            return atoms, alat*Bohr

def read_qe_dynmat(filename, natoms, alat):
    q_point = np.array([0.0, 0.0, 0.0])
    eigval = np.array([])
    eigenvec = np.array([])
    output = open(filename)
    while True:
        line = output.readline()
        if not line:
            break

        p = re.search(r'\s*q\s+=\s+([-\d\.]+)\s+([-\d\.]+)\s+([-\d\.]+)', line)
        if p:
            q_point = np.array([float(p.group(1)), float(p.group(2)), float(p.group(3))])
            print "Found q point: ", q_point
            q_point = q_point/alat
            continue # while True

        p = re.search(r'\s*omega.+?([-\d\.]+)\s*\[cm-1\]', line)
        if p:
            eigval = np.append(eigval, [float(p.group(1))], axis=0) # not efficient
            temp = np.array([])
            for i in range(0, natoms):
                p1 = re.search(r'(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)', output.readline())
                if not p1:
                    break

                temp = np.append(temp, [ float(p1.group(1))+1j*float(p1.group(2)), float(p1.group(3))+1j*float(p1.group(4)), float(p1.group(5))+1j*float(p1.group(6))])

            if len(eigenvec) == 0: # first row of eigenvectors
                eigenvec = temp
            else:
                eigenvec = np.vstack((eigenvec, temp))

    return eigval, np.transpose(eigenvec), q_point










