===================================
Tutorial 13: ASE ONETEP interface
===================================

:Author:  Tom Demeyere
:Date:    August 2023

.. role:: raw-latex(raw)
   :format: latex

This tutorial will guide you through the use of the ASE interface to ONETEP. The Atomic Simulation Environment (ASE) is a set of tools and Python modules for setting up, manipulating, running, visualizing and analyzing atomistic simulations. The ASE interface to ONETEP allows you to set up and run ONETEP calculations from Python scripts.

This tutorial will mainly focus on running ONETEP calculations with ASE, if you are not familiar with ASE as a whole, feel free to consult the mini-tutorial here:

- :download:`learn_ase_in_y_minutes.py <_static/tutorial_13/learn_ase_in_y_minutes.py>`

ASE Configuration File
======================

To run ONETEP with ASE, you need to have a configuration file that specifies the path to the ONETEP binary and the location of the pseudopotentials. There is no default configuration file, so you need to create one.

ASE looks for a configuration file named :code:`config.ini` in the default location :code:`$HOME/.config/ase/`. You can change the location of the configuration file by setting the environment variable :code:`ASE_CONFIG_PATH` to the desired path.

The configuration file should follow the pattern (do not put quotes around the values):

.. code-block:: ini

    [onetep]
    command = mpirun -np 10 -v/path/to/onetep/binary
    pseudo_path = /path/to/pseudos


Replace :code:`/path/to/onetep/binary` with the actual path to your ONETEP binary. If you want to use a different location for the configuration file, you can set the :code:`ASE_CONFIG_PATH` environment variable. For example:

.. code-block:: bash

  $ export ASE_CONFIG_PATH="/path/to/custom/config.ini"

Alternatively, if you don't want to use the configuration file, you can create a :code:`OnetepProfile` directly in your script. For example:

.. code-block:: python

  from ase.calculators.onetep import Onetep, OnetepProfile

  profile = OnetepProfile(
    command="mpirun -np 10 -v /path/to/onetep/binary",
    pseudo_path="/path/to/pseudos"
  )
  calc = Onetep(profile=profile)

This will override the configuration file and use the :code:`OnetepProfile` object instead.

Pseudopotentials
================

If no pseudopotentials are passed ASE will try to guess the files based on the element used and the pseudo_path variable. Otherwise you can pass a dictionary with the element as key and the path to the pseudopotential as value. The path can be either absolute or relative to the :code:`pseudo_path`.

.. code-block:: python

    # Explicitly providing each path:
    calc = Onetep(pseudopotentials = {'H': '/path/to/pseudos/H.usp', 'O': '/path/to/pseudos/O.usp'})
    # Using relative paths:
    calc = Onetep(pseudopotentials = {'H': 'H.usp', 'O': 'O.usp'})
    # ASE will try to guess them if you don't provide them:
    calc = Onetep()

For ASE to correctly guess the pseudopotentials, it is best to use a :code:`pseudo_path` that contains only one pseudopotential file for each element. If there are multiple files for the same element, ASE will not be able to guess which one to use.

ONETEP Calculator
=================

Simple calculations can be setup calling the Onetep calculator without any parameters,
in this case ONETEP's default parameters will be used. For more complex cases using the
:code:`keywords` parameters is necessary. The :code:`'keywords'` parameters is a dictionary, in which each of the keys is a string that should be a ONETEP keyword, and the corresponding value is what you want to set that keyword to in the input.

.. code-block:: python

    from ase.calculators.onetep import Onetep

    # Default parameters
    calc = Onetep()

    # Custom parameters
    keywords = {
        'xc' : 'PBE',
        'do_properties' : True,
        'cutoff_energy' : 35,
        'output_detail': 'verbose',
        'elec_energy_tol': 1.0e-5,
    }

    calc = Onetep(keywords=keywords)

Alternatively you can read an already existing input_file with the function :code:`read_onetep_keywords`

.. code-block:: python

    from ase.io.onetep import read_onetep_keywords

    keywords = read_onetep_keywords('input_file.dat')

    # Let's change one specific keyword
    keywords['xc'] = 'PBE0'

    calc = Onetep(keywords=keywords)

Examples
========

Here is an example python script which sets up a calculation on a water molecule:

.. code-block:: python

    from ase.build import molecule
    from ase.calculators.onetep import Onetep

    water = molecule('H2O', vacuum=10)
    
    calc = Onetep(xc='PBE', paw=True)
    water.calc = calc

    water.get_potential_energy()

Here is a more complex example, this time setting up a :math:`\mathrm{Pt}_{13}` cluster and running a geometry optimisation, note that here as far as ONETEP is concerned we are running singlepoint calculations, the geometry optimisation is done by ASE's BFGS optimiser:

.. code-block:: python

    import numpy as np

    from ase.build import molecule
    from ase.calculators.onetep import Onetep
    from ase.cluster import Octahedron
    from ase.optimize import BFGSLineSearch
    
    # Pt13 from ase.cluster
    nano = Octahedron('Pt', 3, 1)
    nano.center(vacuum=10)

    # ONETEP default are atomic units, one can specify 'cutoff_energy' : '600 eV' if needed.
    keywords = {
        'xc' : 'rpbe',
        'do_properties' : True,
        'cutoff_energy' : 35,
        'output_detail': 'verbose',
        'elec_energy_tol': 1.0e-5/len(atoms),
        'edft': True,
    }

    # append = True will not overwrite file at each step
    calc = Onetep(
        append = True,
        keywords = keywords)

    nanoparticle.calc = calc

    opt = BFGSLineSearch(atoms = nano)
    opt.run(fmax=0.1)


Here is an example of setting up an EELS and LDOS calculation on an N-substituted graphene sheet,
demonstrating several more advanced functionalities (tags, species groups, and overrides to
pseudopotentials and atomic solver strings)

.. code-block:: python

    import numpy as np

    from ase.build import graphene_nanoribbon
    from ase.calculators.onetep import Onetep
    from ase.io import write

    sheet = graphene_nanoribbon(10, 10, type='zigzag', vacuum = 10)

    # Get all distances to center of mass
    com = sheet.get_center_of_mass()
    distances_to_com = np.linalg.norm(sheet.positions - com, axis = 1)

    # Find atoms close to com and change one randomly to N
    p, = np.where(distances_to_com < 5)
    to_nitro = choice(p)
    sheet[to_nitro].symbol = 'N'

    shell_rad = np.array([1.5, 2.5, 3.0, 4.0, 4.5])

    tags = np.zeros(len(sheet), dtype=np.int32)

    # We want to tag atoms that are close to the introduced nitrogen
    for idx, rad in enumerate(reversed(shell_rad)):
        # All distances N-C
        dist = norm(sheet[to_nitro].position - sheet.get_positions(), axis = 1)
        # Which ones are closest to rad?
        p, = np.where(dist < rad)
        # Cannot be the nitrogen itself
        p = p[p != to_nitro]
        # Tags them
        tags[p] = len(shell_rad) - idx

    sheet.set_tags(tags)

    tags = ['' if i == 0 else i for i in tags]

    species = np.unique(np.char.add(sheet.get_chemical_symbols(), tags))

    keywords = {
        'species_core_wf' : ['N /path/to/pseudo/corehole.abinit'],
        'species_solver' : ['N SOLVE conf=1s1 2p4'],
        'pseudo_path': '/Users/tomdm/PseudoPotentials/SSSP_1.2.1',
        'xc' : 'PBE', 
        'paw': True,
        'do_properties': True,
        'cutoff_energy' : '500 eV',
        'species_ldos_groups': species,
        'task' : 'GeometryOptimization'
    }

    calc = Onetep(
        keywords = keywords
    )

    # Checking the input before running the calculation  
    write('to_check.dat', sheet, format='onetep-in', keywords = keywords)
 
    sheet.calc = calc
    # Will actually run the geometry optimisation
    # using ONETEP internal BFGS
    sheet.get_potential_energy()


Quickly restart with solvation effect using the soft sphere solvation model:

.. code-block:: python

    from ase.io import read
    from ase.io.onetep import get_onetep_keywords

    # Read from the previous run...
    optimized_sheet = read("onetep.out")

    # Function to retrieve keywords dict from input file...
    keywords = get_onetep_keywords('onetep.dat')
    
    # We add solvation keywords
    keywords.update(
        {
        'is_implicit_solvent': True, 
        'is_include_apolar': True,
        'is_smeared_ion_rep': True,
        'is_dielectric_model': 'fix_cavity',
        'is_dielectric_function' : 'soft_sphere'
        }
    )

    optimized_sheet.calc = Onetep(keywords=keywords)
    optimized_sheet.get_potential_energy()

Important note
==============

If you are not keen about using ASE to run ONETEP calculations, it is always possible to use ASE to write ONETEP input files and run them manually. This should be done by using the general ASE IO modules :code:`ase.io.write` and :code:`ase.io.read` to write and read ONETEP input files. In every examples above, all you need to do is to replace the :code:`get_potential_energy()` call by a :code:`write` call to write the input file, such as :code:`write('input_file.dat', atoms, format='onetep-in', keywords=keywords)`. You can then run the ONETEP binary manually as you always do.

How to use ASE on HPCs
======================

If the HPC you are using has a module system, you can load the conda module and create an environment with the required packages. If you don't have access to a module system, you can install miniforge in your home directory and create an environment there. A tutorial to do so is available at the end of this document.

How does python launch ONETEP under the hood?
---------------------------------------------

When you run a python script with ASE and ONETEP, ASE will both construst the command to be launched and the input file. The command will be constructed based on the :code:`command` key in the ASE configuration file. Or based on the :code:`command` key in the :code:`OnetepProfile` object if you send the profile manually. The command will be executed with the `subprocess` module using the :code:`check_call` function. The inner working of the :code:`check_call` function is to run the command in a subprocess and wait for it to finish. If the command fails, an exception will be raised. To run the command no new shell is created, and all the environment variables are inherited from the parent process. All stdout and stderr will be redirected to the onetep.out and onetep.err files.

The input file will be created by the IO functions of ASE, namely :code:`ase.io.onetep.write_onetep_input`. This function will write the input file in the format expected by ONETEP. This will be automatically done if a calculation is launched via :code:`atoms.get_potential_energy()` or else.

General case
------------

There are two ways to submit job using ASE on HPC, you can directly sbatch the python script by putting the correct shebang at the top of the script, or you can use an additional bash script to submit the job. The bash script will have to activate the environment and run the python ASE script. Here is an example of such a script:

.. code-block:: bash

  #!/bin/bash
  #SBATCH --job-name=ASE_ONETEP
  ...

  conda activate myenv

  module load ... # Load all the modules needed by ONETEP
  export ... # Set all the environment variables needed by ONETEP

  export ASE_CONFIG_PATH="/path/to/scratch/.ase_config.ini"

  python my_ase_script.py

.. code-block:: python

    # Your python script can look like this

    from ase.build import molecule

    from ase.calculators.onetep import Onetep

    water = molecule('H2O')

    keywords = {
        'xc' : 'PBE',
        'do_properties' : True,
        'cutoff_energy' : 35,
        'output_detail': 'verbose',
        'elec_energy_tol': 1.0e-5/len(water),
    }

    calc = Onetep(keywords=keywords)

    water.calc = calc
    water.get_potential_energy()

Make sure that the ONETEP command being used contain :code:`srun` for example: :code:`command = srun /path/to/onetep/binary`. Otherwise the job will not dispatch correctly on the compute nodes. This is no different from launching a normal job, with the expection that ASE takes care of the input file and the command to be launched.

Archer2
-------

Archer2 is a Cray system, and the conda module is **not** available. You should install it by having a look at the instruction at the end of this document. **One of the Archer2's particularity to keep in mind is that compute nodes only have access to the scratch space and not to the home directory.** You should make sure that every files which will be used during the calculation is accessible from the scratch space, most likely this will be: the input files, the pseudopotentials, the executable and conda. This also means that if you are using the ase config file, you should make sure to change its location with the :code:`ASE_CONFIG_PATH` environment variable to the scratch space. Once this is done you should have a working environment to run ASE on Archer2.

Iridis5
-------

Iridis5 is an Intel based HPC, with conda available as a module. You can alternavely install your own Conda if you want following the instruction at the end of this document if you want it. There is no particularity to keep in mind when running ASE on Iridis5, you can use the conda module to create an environment with the required packages. You can then submit a job with the python script directly or with a bash script as shown above. Make sure to use :code:`srun` in the command to dispatch the job on the compute nodes.

Young
-----

The only particularity of Young is that :code:`srun` is not available, instead a home-made wrapper around `mpirun` is made avaible (`gerun`). **This will not cause limitations as long as you keep each job to serial execution.** For example, if you use the ASE NEB module with threading, i.e. launching multiple ONETEP in parallel in the same PBS job, gerun will most likely not distribute the job correctly, and the calculation will either fail, or be very slow. The only way around this is to make use of the :code:`mpirun` command directly and specifying the node to use for each job. Which will not be detailed here, you should probably use another HPC for this kind of calculation.

Other python packages
=====================

Other packages that can be used with Onetep + ASE are numerous, here we do mini-tutorials for some of them.

DFTD3/DFTD4
-----------

DFTD3 and DFTD4 are dispersion correction methods that can be used with ONETEP. These packages also interface with ASE, which is why they can be used in conjunction with ONETEP. To install DFTD3 or DFTD4, you can use the conda package manager. Here is how to install them:

.. code-block:: bash

  conda install -c conda-forge dftd3-python
  conda install -c conda-forge dftd4-python

If you really care about the performance you should probably compile them yourself, although the performance gain should probably be minimal. After installation they can be used in the ASE calculator as follows:

.. code-block:: python

    from ase.build import molecule
    from ase.calculators.mixing import SumCalculator
    from ase.calculators.onetep import Onetep
    from dftd4.ase import DFTD4
    
    atoms = molecule('H2O')

    calc = SumCalculator([DFTD4(method="PBE"), Onetep(xc="PBE")])
    atoms.calc = calc

    atoms.get_potential_energy()

For DFTD3 the code is pretty much the same, just replace :code:`DFTD4` by :code:`DFTD3`. The DFTD3 version requires to have :code:`method` and :code:`damping` parameters set at all time. With both version you can pass an additional parameter :code:`params_tweaks` where you can manually override the internal D3 parameters, see the documentation for more information.

Alloy Catalysis Automated Toolkit (ACAT)
----------------------------------------

ACAT (https://gitlab.com/asm-dtu/acat) is a python package that can be used to automate the setup of ONETEP calculations for (alloy) catalysis. ACAT can be used in conjunction with ASE, and can be installed using pip:

.. code-block:: bash

  pip install acat

The package allow many operations on both surfaces and nanoclusters, the two main classes are the
:code:`ClusterAdsorptionSites` and the :code:`SlabAdsorptionSites`. Which are used to detect all possible binding sites of your systems. Here is a complete example to create onetep input files for an alloyed nanocluster:

.. code-block:: python

    from pathlib import Path

    from acat.adsorption_sites import ClusterAdsorptionSites
    from acat.build.action import add_adsorbate_to_site
    from ase.cluster import Octahedron
    from ase.io import write

    calc_dir = Path("alloy_project_tutorial")
    calc_dir.mkdir(exist_ok=True)

    atoms = Octahedron("Ni", length=7, cutoff=2)

    # Let's create our alloy
    for atom in atoms:
        if atom.index % 2 == 0:
            atom.symbol = "Pt"

    atoms.center(vacuum=5.0)

    # We create the ACAT object with our parameters,
    # Many more are available, check the documentation
    cas = ClusterAdsorptionSites(
        atoms,
        composition_effect=True,
        label_sites=True,
        surrogate_metal="Ni",
    )

    # Only unique sites, we don't want to duplicate calculations
    sites = cas.get_unique_sites(unique_composition=True)

    for site in sites:
        # add_adsorbate_to_site is modifies the object in place
        # so we copy it to avoid modifying the original object
        tmp = atoms.copy()

        add_adsorbate_to_site(tmp, "O", site)

        # We create a unique custom label based on the information
        label = (f"{tmp.get_chemical_formula(mode='metal').lower()}"
                f"_{site['surface']}_{site['site']}_{site['label']}")

        # The directory for this specific calculation
        current_dir = calc_dir / label
        current_dir.mkdir(exist_ok=True)

        # ASE can of course, write onetep input files
        # In practice you would have to specify keywords and pseudopotentials
        write(current_dir / "onetep.dat", tmp, format="onetep-in")


You will have a directory called `alloy_project_tutorial` with a subdirectory for each adsorption site, each containing an input file for ONETEP. You can then run these input files manually or with ASE as shown in the previous examples. Alternatively you can visualise them using the :code:`ase gui` tool.

Phonopy
-------

Phonopy (https://github.com/phonopy/phonopy) is a python package that can be used to calculate phonon properties of materials. and can be installed using pip or conda:

.. code-block:: bash

  pip install phonopy

Phonopy can be used to calculate the phonon band structure of a material. Usually everything is done using the CLI but I personnaly prefer to use the API directly, here is an example for a water molecule:

.. code-block:: python

    from ase.build import molecule
    from phonopy import Phonopy
    from phonopy.structure.atoms import PhonopyAtoms

    from ase.calculators.onetep import Onetep

    water = molecule('H2O', vacuum=10)

    calc = Onetep()

    phonopy_atoms = PhonopyAtoms(symbols=water.get_chemical_symbols(),
                                 positions=water.get_positions(),
                                 cell=water.get_cell())

    phonopy = Phonopy(phonopy_atoms, supercell_matrix=[[1, 0, 0], [0, 1, 0], [0, 0, 1]])

    phonopy.generate_displacements(distance=0.01)

    displacements = phonopy.supercells_with_displacements

    forces = []

    for i, disp in enumerate(displacements):

        disp_dir = Path(f"displacement_{i}")
        disp_dir.mkdir(exist_ok=True)

        atoms = Atoms(disp.get_chemical_symbols(),
                      disp.get_positions(),
                      cell=disp.get_cell()
        )

        calc.directory = str(disp_dir)
        
        atoms.calc = calc

        forces.append(atoms.get_forces())

    phonopy.forces = forces
    phonopy.produce_force_constants()

    phonon.save("ifc.yaml", settings={'force_constants': True})

    print(phonon.get_frequencies_with_eigenvectors((0, 0, 0))[0]*33.356)

With the annoying fact that the :code:`Atoms`` object have to be manually transfered to :code:`PhonopyAtoms` back and forth. The phonon frequencies are in THz, to convert them to cm-1 you have to multiply by 33.356. The `ifc.yaml` file can be used for further processing. See the phonopy documentation for more information.

Many more
---------

There are many more packages that can be used with ONETEP and ASE. Some of them are listed below:

- **pymatgen**: A python package for materials analysis, which can be used to generate structures, calculate band structures, and much more. (https://github.com/materialsproject/pymatgen)
- **phono3py**: A python package for calculating phonon lifetime and thermal conductivity. (https://github.com/phonopy/phono3py)
- **HiPhive**: A python package to compute higher order force constants without using a specific set of configurations. (https://hiphive.materialsmodeling.org/index.html)
- **Sella**: Sella is a utility primarily intended for refining approximate saddle point geometries. Interfaces well with ASE. (https://github.com/zadorlab/sella)
- **QuAcc**: The Quantum Accelerator (QuAcc) is a python package that can be used to create automated workflows and run them concurrently with workflow managers like Parsl, Dask or Covalent. ONETEP has an interface and a few recipes. (https://github.com/Quantum-Accelerators/quacc)

Conda for the Impatient
=======================

Why Conda?
----------
- **Do not pollute your system-wide python, you might regret it**: Conda creates isolated environments, keeping your system Python clean and preventing conflicts between different projects.
- **Stop compiling your tools use binaries by Conda**: Conda can manage packages for various languages, including R, C++, and Fortran, making it a versatile tool for scientific computing.
- **Complement Conda with pip**: While Conda handles most python package installations, you might occasionally need pip for packages not available in Conda repositories.
- **Conda is self-contained**: Install it everywere, no need for root access. Even HPC systems encourage the use of Conda. Conda will not break your system, and you can remove it easily.

Installing Mambaforge on Linux
------------------------------
1. Download the Mambaforge installer (Linux x86_64) from the Conda Forge repository:
  
  ``wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh``

2. Run the installer:
  
  ``bash Mambaforge-Linux-x86_64.sh``

3. Follow the prompts, agreeing to the license and choosing the installation location.

4. Initialize Mambaforge by running:
  
  ``conda init``

5. Close and reopen your terminal for the changes to take effect.

Installing Conda on Windows
---------------------------

To install Conda on Windows, follow these steps:

1. Visit the official Anaconda website (https://www.anaconda.com) and download the Anaconda Navigator.
2. Run the installer and follow the installation prompts. Make sure to select the option to add Conda to your system's PATH environment variable.
3. Once the installation is complete, open the Anaconda Navigator application to manage packages and environments. You can create environments, install packages, and launch Jupyter notebooks directly from the Navigator interface.
4. If you want to install python's package that are only available through pip you can launch a terminal from the navigator inside the environment you want to install the package and run `pip install package_name`

Creating and Managing Environments
----------------------------------
Create a new environment:
 ``conda create --name myenv``

Activate the environment:
 ``conda activate myenv``

Deactivate the environment:
 ``conda deactivate``

Installing Packages
-------------------
Install packages in the active environment:
 ``conda install numpy pandas``

For packages not available in Conda repositories, use pip:
 ``pip install somepackage``

Updating and Removing Packages
------------------------------
Update a package:
 ``conda update somepackage``

Remove a package:
 ``conda remove somepackage``

Update all packages in the current environment:
 ``conda update --all``

Managing Environments
---------------------
List all environments:
 ``conda env list``

Remove an environment:
 ``conda env remove --name myenv``

Export an environment to a YAML file:
 ``conda env export > environment.yml``

Create an environment from a YAML file:
 ``conda env create -f environment.yml``