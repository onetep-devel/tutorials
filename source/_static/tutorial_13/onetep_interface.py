# Here is how to use ASE to launch simple ONETEP calculations.

from ase.calculators.onetep import Onetep
from ase.build import molecule
from os import environ

keywords = {'output_detail': 'verbose',
            'forces_output_detail': 'verbose',
            # ^^^^^^^^^^^^^^^^^^^ if your version of ONETEP is recent enough
            'write_tightbox_ngwfs': False,
            'write_denskern': False}

# TODO: Change this path to the location of your ONETEP binary
onetep_path = '...'

# TODO: Change this to whatever path you have for the pseudopotentials
pseudo_path = 'Pseudo'

# environ allows to set environment variables like in a shell
# Here we set the number of threads to 1
# and the command to launch ONETEP
environ['OMP_NUM_THREADS'] = '1'
environ['ASE_ONETEP_COMMAND'] = \
    'mpirun -np 3 ' \
    + onetep_path


# TODO: create a water molecule using:
# https://wiki.fysik.dtu.dk/ase/ase/build/build.html#ase.build.molecule
water = molecule('H2O')

# TODO: center the molecule in a cell with a 5.0 Å padding using:
#  https://wiki.fysik.dtu.dk/ase/ase/atoms.html#ase.Atoms.center
water.center(vacuum=5.0)


# Test calculation, we do that with a small ngwf radius
calc = Onetep(label='water',
              pseudo_path=pseudo_path,
              keywords=keywords,
              ngwf_radius=6.0
              )

# TODO: set the calculator of the water molecule to the ONETEP calculator
# https://wiki.fysik.dtu.dk/ase/ase/atoms.html#ase.Atoms.set_calculator
water.calc = calc

# The calculation is launched at this point
print('total energy:', water.get_potential_energy())

# We now want to optimize the geometry of the water molecule
# using the internal optimizer of ONETEP

calc = Onetep(label='water_geom',
              pseudo_path=pseudo_path,
              ngwf_radius=6.0,
              keywords={**keywords,
                        **{'task': 'GeometryOptimization',
                           'geom_force_tol': '1.0e-2 eV/Ang'}
                        }
              )

water.set_calculator(calc)
# Launch the geometry optimisation
forces = water.get_potential_energy()