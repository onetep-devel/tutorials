==================================================
Tutorial 1: Setting up Simple ONETEP Calculations
==================================================

:Author: Nicholas D.M. Hine and Chris-Kriton Skylaris
:Date:   July 2024

.. role:: raw-latex(raw)
   :format: latex

Input files
===========

Setting up a ONETEP job involves creating a main input file with the
suffix :code:`.dat` which contains all the required information to
describe both the system and the parameters of the job. This requires
the user to provide input in the form of keywords and blocks. Keywords
are written in the form

::

    keyword: value [unit]

For example, to specify that the task we wish the code to perform is a
Single-Point energy calculation, we would add:

::

    task : SinglePoint

to our input file (note that capitalisation is irrelevant).

If we wish to specify a cutoff energy of 500 eV for our standard grid,
we would add:

::

    cutoff_energy : 500.0 eV

The value in eV's will be converted internally to atomic units (Eh in
this case). If a keyword is not specified in the input file, it is given
a default value which is intended to work across a broad range of
systems. A full list of keywords and blocks, giving their meaning,
syntax and default values, can be found on the ONETEP keyowrd database:
https://onetep.org/resources/keyword-database/.

Blocks are used to define the values of input parameters which need to
contain multiple records, such as the definition of the unit cell.
They take the form:

::

    %block blockname
    a1 a2 a3
    b1 b2 b3
    ...
    %endblock blockname

Most blocks tend not to have a meaningful default value, and must be
specified if the related functionality is to be used. Comments can be
added to input files using the :code:`# or !` characters. Anything after
these characters on a given line will be ignored.

Setting up the Input File
-------------------------

We will start by running a simple job on a silane molecule
SiH\ :math:`_4`. Create a working directory in which to run ONETEP

::

    > mkdir silane
    > cd silane

Create a new input file called `silane.dat` in your favourite text
editor e.g.

::

    > vi silane.dat &

You might like to put a comment at the top explaining what this input
file is for e.g. 

::

     # Simple ONETEP input file for a silane molecule

The first thing is to specify the simulation cell.
The simplest choice is a cubic box with sides of about 40.0 Bohr.
Enter the 3-component cell vectors, one per line, between the :code:`%block
lattice_cart` and :code:`%endblock lattice_cart` keywords.

Second, the atomic species need to be specified, in this case silicon
and hydrogen.
This information needs to be provided between :code:`%block species` and
:code:`%endblock species` keywords.
In this block, we need to specify five pieces of information per
species, separated by spaces:

1. *Your* symbol for the atomic species (this can be the same as the element symbol).
2. The element symbol itself.
3. The atomic number Z.
4. The number of NGWFs per atom.
5. The NGWF radius.

The number of NGWFs required can usually be judged from the symmetry of
the atomic orbitals involved: In this case four for silicon and one for
hydrogen will be adequate (can you think why? Answer: one s- and three 
p-orbitals for each carbon atom and one s-orbtial for each hydrogen atom).

For this molecule, 6.0 Bohr should be a reasonable starting point for
the NGWF radii.
Each atomic species in our calculation needs a pseudopotential file.
The pseudopotential files are specified between :code:`%block species_pot` and
:code:`%endblock species_pot` keywords. You can use the :code:`hydrogen.recpot`
and :code:`silicon.recpot` files from the ONETEP pseudopotentials :code:`pseudo`
directory. Copy them to your working directory now (or make a symbolic
link by :code:`ln -s /path/to/hydrogen.recpot hydrogen.recpot`).

Next, we need to specify the atomic positions, between :code:`%block
positions_abs` and :code:`%endblock positions_abs` keywords.
There is one line per atom. Remember to use your symbol for the atomic
species as defined in the :code:`species` block. The coordinates are assumed
to be given in Bohr unless specified otherwise. While it is not
requirement in ONETEP that all the atoms should lie within the
simulation cell, it is best (for visualisation purposes) to start by
placing the silicon atom at the centre of the cell.
The Si-H bond length is about 2.76 Bohr and silane is a tetrahedral
molecule. The simplest way to work out the coordinates is to note that
tetrahedral bonds can be chosen to lie along unit vectors
(a, b, 0), (−a, b, 0), (0, −b, a) and (0, -b, −a) where a = 2/3 and
b = 1/3.
For example, the vector for the first Si-H bond is (2.2535, 1.5935,
0.0000) Bohr. Add these offsets to the position of your silicon atom
to create the SiH\ :math:`_4` molecule.

The last essential parameter to specify is the kinetic energy cutoff
parameter for the psinc basis set. A reasonable value to start with is
300 eV. Use the :code:`cutoff_energy` keyword and remember to specify the
energy unit as well as the value.

Running the Job
---------------

Examine the output: if you have followed these instructions it should
converge very quickly (8 iterations) to a total energy of around
-6.1897 Eh.

Convergence Convergence Convergence
-----------------------------------

Just as in any form of traditional DFT, we must ensure that our
calculation results are converged with respect to the size of the
basis. In ONETEP, convergence with basis size is controlled by a small
number of parameters, with respect to which the total energy is
variational. In this context, that means the total energy at a given
value of the parameter will be an upper bound to the true, converged
total energy, and increasing the parameter will monotonically decrease
the total energy, which asymptotically tends to its converged value.

Cutoff Energy
~~~~~~~~~~~~~

The first parameter will be familiar to anyone who has carried out
plane-wave DFT calculations: the cutoff energy. This specifies the
kinetic energy of the maximum G-vector of the reciprocal-space grid,
and therefore the spacing of the real-space grid. With a 40 Bohr cell
and a 300 eV cutoff, ONETEP will have chosen a
48 :math:`\times` 48 :math:`\times` 48 grid,
hence a grid spacing of 0.833 Bohr. This may be too coarse: move your old
output file to a new name (e.g., :code:`SiH4.out_Ec300`) and try changing the
cutoff energy to 500 eV, then re-run the job script. You may wish to
add :code:`output_detail: VERBOSE` to your input file, to see exactly what
grids are being used at each cutoff.

Comparing the two outputs, you should see that the total energy has
decreased by around 0.03 Eh (nearly 1 eV, or 0.2 eV/atom). This suggests
300 eV was too low initially. Try increasing the cutoff in steps of
100 eV (You may wish to automate this, by having a loop in your job
script in which the input file is updated and the job run for each
update, if you are suficiently familiar with bash scripting)

Plot the total energy (ET) as a function of cutoff energy. You should
see a monotonic decrease in ET as a function of Ecut: try to evaluate
at what value you think the total energy is converged to about
0.1 eV/atom of its asymptotic limit. Note that the calculation time
increases rapidly with cutoff energy, because the number of grid
points in each FFT box is growing rapidly with cutoff energy, and thus
each FFT takes longer, so do not try going beyond around 1200 eV.

In few cases in reality do we require strict convergence of the total
energy. It is more usual that we require convergence of some
measurable quantity such as a binding energy, which is based on energy
differences.
In that case, we do not require the total energy to be converged, only
the difference between total energies of very similar systems. This
may converge much faster than the total energy itself, presuming the
same species are present in both systems. Always consider what it is
that you need converging before you start running enormous
calculations!

NGWF radius
~~~~~~~~~~~

Next, we will investigate convergence with respect to the NGWF radius.
Pick a value of cutoff energy for which you can perform reasonably
fast calculations (say, 500.0 eV) and try increasing the NGWF radius
from 6.0 to 10.0 in 1.0 Bohr steps. Plot the total energy against NGWF
radius. Again, you should see a monotonic decrease. Note that above
6.0 Bohr the FFT box is as large as the simulation cell, in a larger
cell this would keep growing, and the calculation time would increase
rapidly. Also, you should notice that the number of NGWF Conjugate
Gradients iterations grows with the size of the localisation region,
this is natural since with larger spheres there are more NGWF
coeficients to simultaneously optimise.
You may also wish to try converging with respect to the number of
NGWFs per atom (e.g., try 9 NGWFs on the Silicon). In some systems,
notably crystalline solids, this can be crucial to achieving good
convergence of the NGWFs themselves.

Kernel Cutoff
~~~~~~~~~~~~~

This SiH\ :math:`_4` system is too small to investigate convergence with
respect to the cutoff of the density kernel. In larger systems
truncation of the density kernel can be a good way to speed up the
calculation. Indeed, asymptotically it is only by truncating the kernel
that true 'linear-scaling' behaviour of the computational effort will
be observed.

The kernel cutoff is controlled by the :code:`kernel_cutoff` keyword. This
defaults to 1000 Bohr (i.e. effectively infinite). Density kernel
truncation should be used with a degree of caution: generally speaking,
one would want to be able to run a full calculation for a fairly large
system first, with an infinite cutoff, to establish a known baseline.
Then, try decreasing the kernel cutoff from that point and see what the
effect is on the total energy, on the level of NGWF convergence (as
measured by the NGWF RMS gradient), and on the computation time. If
significant time savings can be achieved without trading in too much
accuracy, it may be worthwhile to bring down the cutoff for all similar
calculations. Proceed with care, though as calculations with a truncated
kernel tend to converge in a less stable manner.

Crystalline Silicon
-------------------

You may wish to try out also a calculation on a periodic solid.
As it is fairly well-behaved but illustrates some interesting
concepts, let's try crystalline silicon, in the diamond (f.c.c.)
structure. We will build a 2 :math:`\times` 2 :math:`\times` 2 version of the
8-atom simple-cubic unit cell.
Copy your SiH\ :math:`_4` input to a new file (e.g., :code:`Si64.dat`) in a new
directory (e.g., :code:`SILICON`) and remove the references to hydrogen from the
species and :code:`species_pot` blocks. Copy :code:`silicon.recpot` to this
directory as well. In the new input file, set the NGWF radius to 7.0 Bohr, the
number of NGWFs per atom to 4, and the cutoff energy to 600 eV. Edit
the cell side length so that it is 2 :math:`\times` the lattice parameter of
crystalline silicon in the LDA (around 10.1667 Bohr). For reasons that
will become clear if you read the last section of this tutorial, on
Common Problems, also set :code:`ngwf_cg_max_step: 8.0` to prevent the CG line
step being capped unnecessarily and maxit_ngwf_cg: 30 to terminate the
NGWF CG after 30 iterations (in case it's not converging).
Typing out the positions would be rather time-consuming and
error-prone with 64 atoms in the cell, so use your favourite
scripting/programming language (bash, awk, python, perl, etc would all
be suitable or even C or FORTRAN) to write a list of the positions.
You will need to repeat the basis (atoms at (0, 0, 0) and
(1/4,1/4,1/4)a ) at each of the positions of the f.c.c. lattice: (0,
0, 0), (1/2,1/2,0), (0,1/2,1/2), and (1/2, 0, 1/2,). Copy the result
into your :code:`positions_abs` block. An example input file for this job can
be found on the tutorial web page.
The calculation should take around a few minutes. Feel free to stop it
as soon as you see what is happening, since you will find that the
calculation fails to converge: the RMS gradient remains stuck above
the threshold for convergence. Likewise, the total energy will not
converge to a fixed value. Make a copy of your output and modify the
NGWF radius in the input file to 8.0 Bohr and the number of NGWFs per
Si atom to 9. This introduces NGWFs with d-like symmetry rather than
just s and p, allowing much more variational freedom. You should now
find the calculation converges nicely, but will take rather longer to
run.

Now try activating :code:`write_forces: T` to calculate the forces on each
atom. All the forces should be small: in principle they are
constrained by the symmetry of the crystal to be exactly zero.
However, you will see that they are not exactly zero because the
symmetry of the system is broken by the psinc grid, which is not
necessarily commensurate with the unit cell. However, in this small
cell, it will not be possible to fix this as the number of points
across the FFT box must be odd, and in a small cell the simulation
cell and the FFT box coincide, so the number of points across the
simulation cell must also be odd.

Adjust your script to write a 5 :math:`\times` 5 :math:`\times` 5 supercell of
the crystal (1000 atoms). Reduce the kernel cutoff to 25 Bohr with
:code:`kernel_cutoff: 25.0` and set the code to perform 1 NGWF iteration only
:code:`maxit_ngwf_cg: 1`
(otherwise the calculation would take longer to run, you can try this if
you have time). To restore the symmetry, adjust the psinc_spacing value
to be a divisor of the supercell length such that an exact number of
grid points spans each unit cell of the crystal (pick a value which
gives an effective cutoff energy close to 600.0 eV so as not to increase
the run time too much) and and set off the 1000 atom job. This should
not take too long on 32 cores.

Beyond around 500 atoms, the calculation should be into the so-called
'linear-scaling' regime, so the 8000 atom calculation should only take
a little over 8 times the 1000 atom calculation. This is rather better
than the nearly 512 times longer it would take with traditional
cubic-scaling DFT!

Diagnosing Common Failures
--------------------------

With badly-chosen input settings, even fairly standard calculations in
ONETEP will not converge, or may even converge to the wrong result.
Fortunately, many of these problems are easy to fix with a bit of
experience. In general, it is advisable to run with full output
verbosity (:code:`output_detail: VERBOSE`) the first few times you run a new
kind of system, and to be on the lookout for any warnings or garbage
numbers in the output (e.g., :code:`****`'s in place of what should be real
numbers). Remember that for the energy to be accurate, we must have
simultaneous convergence of both the density kernel and the NGWFs. If
either of these are not converging well by the end of the calculation,
there may be a problem. In this section, we will briefly examine some
reasons behind common types of convergence failure, and what to do to
eliminate those failures and perform accurate simulations.

- **Problem**: Repeated 'safe' steps (of 0.150 or 0.100) during NGWF
  Conjugate Gradients optimisation, leading to poor or no convergence.
  This often means that the step length cap for NGWF CG is too short.

  | **Solution:** increase :code:`ngwf_cg_max_step`, e.g., to 8.0.

- **Problem**: Repeated 'safe' steps (of 0.150 or 0.100) during LNV
  Conjugate Gradients optimisation, leading to poor or no convergence.
  This often means that the step length cap for LNV CG is too short.

  | **Solution:** increase :code:`lnv_cg_max_step`, e.g., to 8.0.

- **Problem**: Occupancies 'break' during LNV optimisation of kernel.
  Examine the output with :code:`output_detail: VERBOSE` and look at the
  occupancy error and occupancy bounds during the "Penalty functional
  idempotency correction" section of each LNV step.
  Check for occupancies outside the stable range (approx -0.3:1.3) or RMS
  occupancy errors not decreasing (particularly if no kernel truncation is
  applied).

  | **Solution:** Activate
    LNV line step checking with :code:`lnv_check_trial_steps: T`. This checks
    that the kernel is still stable after the proposed line step is taken.

- **Problem**: Occupancies are 'broken' from start of calculation. Symptoms
  as above. Palser Manolopoulos may be unstable due to degeneracy or
  near-degeneracy at the Fermi level. Check the output of Palser
  Manolopoulos for warnings.

  | **Solution:** If there is an initial degeneracy
    at the Fermi level, an O(N3) diagonalisation may be required to get a
    good starting kernel. Set :code:`maxit_palser_mano : -1`.

- **Problem**: RMS Commutator (HKS-SKH) of kernel and Hamiltonian stagnates
  (stops going down with each iteration) during LNV optimisation. This is
  a sign that the current set of NGWFs is not able to represent a density
  matrix that both reproduces the electron density that generated the
  Hamiltonian while simultaneously describing the occupied eigenstates of
  that Hamiltonian. If this problem does not start to go away after a few
  steps of NGWF optimisation, a better or larger initial set of NGWFs may
  be required.

  | **Solutions:** Increase number of NGWFs per atom, increase radius of NGWFs,
    improve initial guess for NGWFs.

- **Problem**: RMS NGWF gradient stagnates (stops going down) during NGWF CG
  optimisation, while energy is still going down slowly. This often
  suggests that the NGWFs may have expanded away from their centres to
  have significant value near the edge of their localisation region, and
  thus cannot optimise successfully.

  | **Solution:** Increase NGWF radius.
    Sometimes increasing the kinetic energy cutoff helps as well. For
    smaller systems and initial tests, a useful check on the accuracy of the
    final result is to perform a full O(N\ :math:`^3`) diagonalisation at
    the end of the calculation, if it is computationally feasible to do so.
    To activate this, turn on a properties calculation with :code:`do_properties: T`
    , and then ask for an eigenvalue calculation of the first 100
    eigenvalues either side of the Fermi energy, for the kernel and
    Hamiltonian matrices, by setting :code:`num_eigenvalues: 100`. If all is
    well, then the occupation eigenvalues should all be close to 0.00000 or
    1.00000 (empty or full) and the Hamiltonian eigenvalues should all be
    within a sensible range.

One final note if you're not getting the result you expect - check the
units on your atomic positions! ONETEP expects positions in Bohr if the
units are not specified, so if your positions are in Angstroms, you will
need to add 'ang' as the first line of the :code:`positions_abs` block.

This completes tutorial 1.

Files for this tutorial:

 - :download:`SiH4.dat <_static/tutorial_1/SiH4.dat>`
 - :download:`Si8000.dat <_static/tutorial_1/Si8000.dat>`
 - :download:`Si64.dat <_static/tutorial_1/Si64.dat>`
 - :download:`Si1000.dat <_static/tutorial_1/Si1000.dat>`


