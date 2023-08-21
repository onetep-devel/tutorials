# Delete this line when you understand the code
#raise NotImplementedError('>>>>>> Please read this file before running it! <<<<<<')

# The Atomic Simulation Environment (ASE) is a set of tools and Python modules for setting up, manipulating, running, visualizing and analyzing atomistic simulations. The code is freely available under the GNU LGPL license.
# https://wiki.fysik.dtu.dk/ase/

#  One of the main object in ASE is the Atom object. It is a container for a single atom and its properties. The properties are stored as attributes of the Atom object. For example, the position of the atom is stored as the attribute position. The position is a three-dimensional vector, and can be accessed as a NumPy array:

from ase import Atom

a = Atom('C', position=(0, 0, 0))

print(a.position)

# Change the position of the atom.
a.position = (1, 2, 3)

# The Atom object is a subclass of the Atoms object. The Atoms object is a container for a collection of Atom. The Atoms object also has a number of attributes, such as the unit cell, and methods for accessing and modifying the attributes.

from ase import Atoms
a = Atoms('CH3CH3CHOH', positions=[(0, 0, 0), (0, 0, 1), (0, 0, 2),
                                 (0, 1, 0), (0, 1, 1), (0, 1, 2),
                                 (1, 1, 0), (1, 1, 1), (1, 1, 2),
                                 (1, 0, 0), (1, 0, 1), (1, 0, 2)])

# ASE also has a set of 'getter' and 'setter' methods for accessing and modifying the attributes. For example, the position can be accessed and modified using the get_position and set_position methods:

# This is equivalent to a.positions
print(a.get_positions())

# It might be useful to center the molecule

a.center(vacuum=10.0)

# Access the cell
print(a.get_cell())
# Access the cell vectors by accessing the cell attribute
print(a.cell.cellpar())
# Access the cell inverse
print(a.cell.reciprocal())
# Access the cell volume
print(a.get_volume())
# Access the center of mass
print(a.get_center_of_mass())
# Access the moment of inertia
print(a.get_moments_of_inertia())
# Access the number of atoms of a given type
print(a.get_atomic_numbers())
# The full list is available here: https://wiki.fysik.dtu.dk/ase/ase/atoms.html#ase.Atoms

# The Atoms object can be sliced and indexed like a list. For example, the first atom can be accessed using a[0], and the first three atoms can be accessed using a[:3]:

print(a[0])
print(a[:3])

# The advantage of using a.positions over the getter and setter methods is that the positions can be accessed and modified in a single operation:

# Change the positions of the first atom or the first three atoms
a.positions[0] = (1, 2, 3)
a.positions[:3] = [(1, 2, 3), (4, 5, 6), (7, 8, 9)]

# It is possible to use more complex operations, any operation that can be applied to a NumPy array can be applied to the positions:

# Translate all atoms by that much
a.positions += (1, 2, 3)

# Calculate the norm of the distance between the first atom and all other atoms (We use numpy broadcasting here)
print(((a.positions[0] - a.positions)**2).sum(axis=1)**0.5)
# Does the same
import numpy as np
print(np.linalg.norm(a.positions[0] - a.positions, axis=1))


####################################################
# 1. Creating Atoms objects
####################################################

# Instead of manually creating atoms object, you can access of bunch of methods.

from ase.build import molecule

water = molecule('H2O')

# Internally ASE uses a database, it will not work for everything:

# propanol = molecule('CH3CH2CH2OH') < -- This will not work, If really needed
# there is a pubchem interface that can be used to download the structure from pubchem
# https://wiki.fysik.dtu.dk/ase/ase/data.html#ase.data.pubchem.pubchem_atoms_conformer_search

# Create a bulk crystal of aluminum with the face-centered cubic (fcc) structure:

from ase.build import bulk

aluminum = bulk('Al')

# Create a slab of aluminum with the fcc structure and the (111) surface:

from ase.build import fcc111

aluminum = fcc111('Al', size=(2, 2, 3), vacuum=10.0)

# Create a graphene sheet:

from ase.build import graphene_nanoribbon

graphene_sheet = graphene_nanoribbon(5, 5, type='armchair', saturated=True, vacuum=3.5)

# Most of the methods for creating atoms objects are listed here: https://wiki.fysik.dtu.dk/ase/ase/build/build.html
# Some other methods are scattered accross the documentation: 
#https://wiki.fysik.dtu.dk/ase/ase/spacegroup/spacegroup.html
#https://wiki.fysik.dtu.dk/ase/ase/cluster/cluster.html

# Because the documentation can be 'scattered', the best way to know if ASE is capable of creating a particular structure is to use google.

####################################################
# 2. Reading and writing Atoms objects
####################################################

# ASE can read and write atoms objects in 50+ formats, the full list is available here: https://wiki.fysik.dtu.dk/ase/ase/io/io.html

# Only two functions are necessary to read and write
# ASE internally call the right function depending on the extension of the file
# Or looking at the content of the file. It might fail sometimes, so you will need to specify the format.
# by passing the 'format' keyword.

from ase.io import read, write

# If you possess an atom object, you can write it to a file:

water.write('water.xyz')

# And read it again...
atoms = read('water.xyz')

# We can write it in different formats...
# Write highly specialized LAMMPS data file:
write('water.lammps-data', atoms, units='real', atom_style='charge', specorder=['O', 'H'], masses=True)

# Write a onetep input file with default settings:
# Might not work if you don't have the new interface.
write('water.dat', atoms, format='onetep-in')

# Want a quick picture?
write('water.png', atoms)

# ASE introduced a new 'extxyz' which extend the xyz format to contains info about the unit cell, energies, forces or other informations. This will be done atomatically if you save to a file with the .xyz extension.

####################################################
# 3. Calculators
####################################################

# ASE has a large number of calculators, the full list is available here: https://wiki.fysik.dtu.dk/ase/ase/calculators/calculators.html

# The simpliest calculator is the lennard-jones calculator, it is a simple pair potential calculator:

from ase.calculators.lj import LennardJones

atoms = molecule('H2O')

atoms.calc = LennardJones()

print(atoms.get_potential_energy())
print(atoms.get_forces())

# The calculator can be modified using the set method:

atoms.calc.set(sigma=3.0, epsilon=0.1)

print(atoms.get_potential_energy())

# Each calculator has a set of parameters that can be modified, the full list is available here: https://wiki.fysik.dtu.dk/ase/ase/calculators/calculators.html#ase.calculators.calculator.Calculator

# Not all calculators are created equal, some of them are not well coded/maintened...

# Has you might have guessed, this can be used to perform high throughput workflows. For example, we can calculate the energy of a large number of structures:

aluminum = bulk('Al', cubic=True)

aluminum.calc = LennardJones()

original_cell = aluminum.cell.copy()
# Monte Carlo strain on the cell
for i in range(10):
    strain = np.random.normal(1, 0.05, 3)
    print('Strain:', strain)
    strain = np.diag(strain)
    aluminum.set_cell(original_cell @ strain, scale_atoms=True)
    print(aluminum.get_potential_energy())

# This is a very simple example, but it can be used to perform more complex workflows.

# You can use phonopy to perform phonon calculations: https://phonopy.github.io/phonopy/phonopy-module.html
# You can use acat to build complex catalytic systems: https://asm-dtu.gitlab.io/acat/build.html
# You can use wulffpack to automatically construct nanoparticles: https://wulffpack.materialsmodeling.org
# Etc... Recently I helped a collegue to code a fully functional force-bias monte-carlo and plug it into ASE, it took less than 100 lines of code.

# You can create all kind of workflow using numpy, scipy, scikit to perform basic or complex operations on atoms objects and the results from the calculators.

####################################################
# 4. Integrated workflows
####################################################

# ASE has a set of integrated workflows, geometry optimisation, molecular dynamics, nudged elastic band, etc...
#https://wiki.fysik.dtu.dk/ase/ase/optimize.html
#https://wiki.fysik.dtu.dk/ase/ase/md.html
#https://wiki.fysik.dtu.dk/ase/ase/neb.html

from ase.calculators.emt import EMT
from ase.md.langevin import Langevin
from ase.units import fs

# Need initial velocities?
# aluminum.set_velocities(np.random.normal(0, 1e-18, (len(aluminum), 3)))

pt_bulk = bulk('Pt', cubic=True)
pt_bulk.calc = EMT()

dyn = Langevin(pt_bulk, timestep=2.0*fs, friction=1e-2, temperature_K=1000.0,
                     trajectory='md.traj', logfile='md.log')

dyn.run(1000)

# In the next section we will visualize the trajectory

####################################################
# 5. Visualizing Atoms objects
####################################################

# ASE has a built-in GUI, it is not the best GUI in the world, but it is good enough for most cases.

# The GUI can be started using the view function:
from ase.visualize import view

view(atoms)

# Many options are possible within the gui, move, rotate, change the cell... have a look at
# the options.

# It is also possible to visualize a file directly from the command line:

# ase gui water.xyz

####################################################
# 6. Tips and tricks
####################################################


# By default the read function will return the last atoms object in the file, if you want to read all the atoms objects in the file, you can use the index keyword:

atoms = read('md.traj', index=':')
atoms = atoms[0]

# Make use of Numpy functions.

# Translate only oxygens
p = atoms.numbers == 8
atoms.positions[p] += 0.01


# Check if two atoms objects are symmetrically the same
from ase.utils.structure_comparator import SymmetryEquivalenceCheck
comp = SymmetryEquivalenceCheck()

a = bulk('Al')
b = a.copy()
b.rotate(60, 'x', rotate_cell=True)
print(comp.compare(a, b))

# Easily create supercells, np.eye(3) returns the identity matrix of size 3x3
# np.eye(3)*5 create a 3x3 matrix with 5 on the diagonal.
from ase.build import make_supercell
c = make_supercell(a, np.eye(3)*5)

# Why does this return False?
print(comp.compare(a, c))

# Because the atoms are not in the same primitive cell, we can fix this by setting to_primitive=True
comp = SymmetryEquivalenceCheck(to_primitive=True)

print(comp.compare(a, c))

# To be continued...