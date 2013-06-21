# striped interface for QE
import re
import fileinput
import numpy as np
from modulation.units import Bohr

def read_qe(filename, symbols=None):
    alat = 0.0

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
            cell = []
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

            print positions
            print expanded_symbols
            #atoms = Atoms(symbols=expaned_symbols,
            #          cell=cell,
            #          scaled_positions=positions)



def get_atoms_from_output(f, symbols):
    lines = f.readlines()

    line1 = [x for x in lines[0].split()]
    if is_exist_symbols(line1):
        symbols = line1

    scale = float(lines[1])

    cell = []
    for i in range(2, 5):
        cell.append([float(x) for x in lines[i].split()[:3]])
    cell = np.array(cell) * scale

    try:
        num_atoms = np.array([int(x) for x in lines[5].split()])
        line_at = 6
    except ValueError:
        symbols = [x for x in lines[5].split()]
        num_atoms = np.array([int(x) for x in lines[6].split()])
        line_at = 7
    
    expaned_symbols = expand_symbols(num_atoms, symbols)

    if lines[line_at][0].lower() == 's':
        line_at += 1

    is_scaled = True
    if (lines[line_at][0].lower() == 'c' or
        lines[line_at][0].lower() == 'k'):
        is_scaled = False

    line_at += 1

    positions = []
    for i in range(line_at, line_at + num_atoms.sum()):
        positions.append([float(x) for x in lines[i].split()[:3]])

    if is_scaled:
        atoms = Atoms(symbols=expaned_symbols,
                      cell=cell,
                      scaled_positions=positions)
    else:
        atoms = Atoms(symbols=expaned_symbols,
                      cell=cell,
                      positions=positions)
        
    return atoms
