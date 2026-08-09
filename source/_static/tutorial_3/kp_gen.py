#!/usr/bin/env python
"""kp_gen.py: Generate a KPOINTS file for ONETEP band structure calculations.

This script generates a list of k-points along the high-symmetry path of the
Brillouin zone for the given crystal structure, using seekpath to determine
the path automatically. The resulting file (KPOINTS.k_path by default) is
meant to be appended to an ONETEP .dat input file to perform a "fake-kpoints"
band structure calculation, in which the k-points along the path are given
zero weight so that they do not affect the self-consistent solution.

The crystal structure may be supplied in any format supported by ASE (e.g.
POSCAR or CIF), or directly as an ONETEP input file, which is detected by the
presence of a %block construct and parsed with ase.io.onetep.read_onetep_in.

Usage
-----
    python kp_gen.py [-h] [-c INPUT_FILE] [-o OUTPUT_FILE] [-r RESOLUTION]
                     [-t] [-s SYMPREC] [-e] [--hybrid] [--vdir VDIR] [-v]

Options
-------
-c, --input <file>       Crystal structure file (POSCAR, CIF, ... or an
                         ONETEP .dat input). Default: POSCAR.
-o, --output <file>      File to which the k-points are written.
                         Default: KPOINTS.k_path.
-r, --resolution <val>   Reference target distance between neighbouring
                         k-points along the path, in units of 1/Ang.
                         Default: 0.1.
-t                       Turn off time reversal symmetry.
-s, --symprec <val>      Precision for symmetry detection (spglib).
                         Default: 0.01 Ang.
-e, --explicit           Write an explicit list of k-points (default). Use
                         -e to switch to line-mode output instead.
--hybrid                 Prepare the k-point list for hybrid band structure
                         calculations (k-point weights set to zero).
--vdir <0|1|2>           Vacuum direction for two-dimensional materials
                         (0/1/2 = x/y/z); segments that move along this
                         direction are removed from the path.
-v, --verbose            Verbose output, including symmetry information and
                         the generated k-point path.

Output
------
By default an explicit list of k-points is written in a format compatible
with the ONETEP %block kpoints_list construct:

    x y z weight !label

If -e is given, the output is written in line-mode with 30 points per
segment instead.
"""

from __future__ import print_function
import seekpath
import sys
import re
import ase
import numpy as np
import ase.io as io
import argparse
import time
from ase.io.onetep import read_onetep_in


parser = argparse.ArgumentParser(description='A script to make KPOINTS file.')
parser.add_argument("-c",
                  action="store", dest="input_file", default="POSCAR",
                  help="The crystal structure. Default: POSCAR")
parser.add_argument("-o", "--output",
                  action="store", dest="output_file", default="KPOINTS.k_path",
                  help="The file to which the new kpoints will be appended. Default: KPOINTS.kpgen.")
parser.add_argument("-r", "--resolution",
                  action="store", dest="resolution", default=0.1,
                  help="a reference target distance between neighboring k-points in the path, in units of 1/ang.")
parser.add_argument("-t", action="store_true", dest="time_reversal",
                  help="Turns off time reversal symmetry.")
parser.add_argument("-s", "--symprec", action="store", default=0.01, dest="symprec",
                  help="precision for symmetry detection [spglib]. Default: 0.01 \AA")
parser.add_argument("-e", "--explicit", action="store_false", default=True, dest="explicit",
                  help="write explicit kpoints? Default: True")
parser.add_argument("--hybrid", action="store_true", dest="hybrid",
                  help="For hybrid bandstructure calculation?")
parser.add_argument("--vdir", action="store", dest="vdir", default=None,
                  help="vacuum dir? [0->x;1->y;2->z]")
parser.add_argument("-v", action="store_true", dest="verbose", default=False,
                  help="verbose output?")
options = parser.parse_args()

# Starting
#----------------------------
if options.verbose == True:
    starttime = time.time()
    print("Starting calculation at", end='')
    print(time.strftime("%H:%M:%S on %a %d %b %Y\n"))

# read in structure
def read_structure(path):
    """Read a crystal structure, supporting both standard ASE formats
    (POSCAR, CIF, ...) and ONETEP input files.

    ONETEP inputs are detected by their '%block'/'%BLOCK' syntax, which
    no other ASE format uses. They are parsed with
    ase.io.onetep.read_onetep_in, which returns a dict whose 'atoms' key
    holds the ase.Atoms object (cell and positions).
    """
    with open(path) as fd:
        head = fd.read(4096)
    if re.search(r"%block", head, re.IGNORECASE):
        if options.verbose == True:
            print("ONETEP input file detected: %s" % path)
        with open(path) as fd:
            atoms = read_onetep_in(fd)["atoms"]
        return atoms
    return io.read(path)


structure = read_structure(options.input_file)
numbers = structure.get_atomic_numbers()
inp = (structure.cell,structure.get_scaled_positions(),numbers)


# Turn off time reversal symmetry if necessary
if not options.time_reversal:
    tr = True
else:
    tr = False

# get K-points
explicit_data = seekpath.get_explicit_k_path(inp,with_time_reversal=tr,
                                             reference_distance=float(options.resolution),
                                             symprec=float(options.symprec))
kpath = explicit_data['explicit_kpoints_rel']
seg = np.array(explicit_data['explicit_segments'])
labels = np.array(explicit_data['explicit_kpoints_labels'])

# return symmetry
if options.verbose == True:
    print("Structure information:")
    print("\tPrecision for finding symmetry: %8.6f \AA" % options.symprec)
    print("\tSpace group number: %s" % explicit_data['spacegroup_number'])
    print("\tSpace international symbol: %s" % explicit_data['spacegroup_international'])
    print("\nTime reversal symmtery: %s" % tr)


# 2D material?
seg_rm = []
if options.vdir!=None:
    for iseg in range(len(seg)):
        if kpath[seg[iseg,0],int(options.vdir)]!=0.0 or kpath[seg[iseg,1]-1,int(options.vdir)]!=0.0:
            seg_rm.append(iseg)
seg = np.delete(seg, seg_rm, axis=0)

# construct path
fkpath = np.array([]).reshape(0,3)
for iseg in range(len(seg)):
    fkpath=np.append(fkpath,kpath[seg[iseg,0]:seg[iseg,1]],axis=0)

# construct label path
flabels = np.array([])
for iseg in range(len(seg)):
    flabels=np.append(flabels,labels[seg[iseg,0]:seg[iseg,1]],axis=0)

if options.verbose == True:
    print("k-point path:")
    for iseg in range(len(seg)):
        print("\t%s\t(%8.6f %8.6f %8.6f)\t->\t%s\t(%8.6f %8.6f %8.6f)" %(labels[seg[iseg,0]],
                                                                         kpath[seg[iseg,0],0],
                                                                         kpath[seg[iseg,0],1],
                                                                         kpath[seg[iseg,0],2],
                                                                         labels[seg[iseg,1]-1],
                                                                         kpath[seg[iseg,1]-1,0],
                                                                         kpath[seg[iseg,1]-1,1],
                                                                         kpath[seg[iseg,1]-1,2]))

# construct label path:
#------------------------
label_path = []
for iseg in range(len(seg)):
    label_path.append(f"""{kpath[seg[iseg,0],0]:+8.6f} {kpath[seg[iseg,0],1]:+8.6f} {kpath[seg[iseg,0],2]:+8.6f} !{labels[seg[iseg,0]]}
{kpath[seg[iseg,1]-1,0]:+8.6f} {kpath[seg[iseg,1]-1,1]:+8.6f} {kpath[seg[iseg,1]-1,2]:+8.6f} !{labels[seg[iseg,1]-1]}\n\n""")

#print(''.join(label_path))

# set k-point weight to zero?
if options.hybrid:
    weight = 0.
    if options.verbose == True:
        print("For hybrid calculations.")
else:
    weight = 1./len(fkpath)

# write data
with open(options.output_file,'w') as outfile:
    if options.explicit:
        outfile.write("File generated by kp_gen.py\n")
        outfile.write(str(len(fkpath))+"\n")
        outfile.write("Reciprocal\n")
        for ipoint in range(len(fkpath)):
            outfile.write("% 8.6f % 8.6f % 8.6f %5.3f !%s\n" % (fkpath[ipoint,0],
                                                             fkpath[ipoint,1],
                                                             fkpath[ipoint,2],
                                                             weight,
                                                             flabels[ipoint]))
    else:
        outfile.write("File generated by kp_gen.py\n")
        outfile.write(str(30)+"\n")
        outfile.write("Line-mode\n")
        outfile.write("rec\n")
        outfile.write(''.join(label_path))

if options.verbose == True:
    print('Output written to {}'.format(options.output_file))
    endtime = time.time()
    runtime = endtime-starttime
    print("\nEnd of calculation.")
    print("Program was running for %.2f seconds." % runtime)
