========================================================
Tutorial 14: Electron Energy Loss Spectroscopy in ONETEP
========================================================

:Author: Edward Tait and Nicholas Hine
:Date:   August 2023 (original text c2018)

.. role:: raw-latex(raw)
   :format: latex
..

Introduction
============

EELS (Electron Energy Loss Spectroscopy) is a spectroscopic technique which 
combines high spatial resolution with fair energy resolution, and is thus
sensitive to local electronic structure in a material. The theory behind
EELS calculations in ONETEP is explained in the 
`documentation page on EELS <https://docs.onetep.org/eels_in_onetep.html>`,
as well as in a paper giving an overview of the capabilities [1], and in the
thesis of one of the authors (Edward Tait) [2]. 

EELS calculations in ONETEP proceed by running first a singlepoint calculation
(with or without a core hole, as required) and then running a properties
calculation. The properties calculation will output files suitable for use
with the OptaDoS package, after a little renaming and rearrangement which can
be carried out with a script provided. 

System setup
============

We will demonstrate the procedure for running an EELS calculation
on a toy system: silene (the silicon equivalent of ethene). An example input
file is provided in the Files section below, as are the PAW potentials for
Silicon and Hydrogen and the associated core wavefunction data for Silicon. 
The tutorial files use the JTH pseudopotentials[5], with the addition
of core orbitals, and we also regenerated the PAW potential and Core
orbitals with a core hole, for the later part of this tutorial.

There are a few differences between the ONETEP input here and one for a
typical single point calculation:

- PAW is mandatory
- We specify a second species for the atom whose core electrons we are exciting
- Because a conduction calculation is being performed we must provide a ``species_cond`` block
- We must provide a ``species_core_wf`` block and every species must be listed there.

The input file can be run swiftly on a single node and should produce a
large number of output files. Most of these files are ``.cube`` files of wavefunctions
produced by default during the properties calculations. The files of interest
are the ``.elnes_bin`` files, which contain OptaDoS compatible matrix elements.

A little more setup is needed before we can run OptaDoS (using the silene
example):

- A dummy castep ``silene-out.cell`` file must be produced, and it must contain a symmetry block
- By default two ``.elnes bin`` files are produced, one based on Kohn-Sham wavefunctions represented using only the valence NGWFs (``silene_val_grad.elnes_bin``) and a second which makes use of the joint basis of valence and conduction NGWFs (``silene_joint_grad.elnes_bin``)
- As per the discussion above, you should choose the latter and copy it to ``silene.elnes bin``.
- A ``silene.bands`` file must be produced, this is best done by copying ``silene.joint bands`` to ``silene.bands``.
- An OptaDoS input file, ``silene.odi`` is needed.

To assist in these tasks a utility script, ``prep_optados_eels``, is provided in the
utils folder of the onetep distribution. Run it with the calculation seed name
as its argument and the steps listed above will be completed automatically.

The .odi file produced should be regarded as a basic template, consult the
OptaDoS documentation[6] if you wish to use more advanced features. Note
that at the moment only fixed broadening is supported by onetep.
When you are satisfied with your OptaDoS input file, execute OptaDoS
with your calculation seed name as the argument. All being well, you should
see a .dat file which you can plot with your favorite tool. Individual edges are
listed sequentially in the file, so a little post processing with awk or python is
needed to separate the edges for individual plotting


Files for this tutorial
=======================

 - :download:`Si2H4_EELS_Example.dat <_static/tutorial_14/Si2H4_EELS_Example.dat>`
 - :download:`H.PBE-paw.abinit <_static/tutorial_14/H.PBE-paw.abinit>`
 - :download:`Si.PBE-paw.abinit <_static/tutorial_14/Si.PBE-paw.abinit>`
 - :download:`Si.PBE-paw.corewf.abinit <_static/tutorial_14/Si.PBE-paw.corewf.abinit>`
 - :download:`Si_corehole.PBE-paw.abinit <_static/tutorial_14/Si_corehole.PBE-paw.abinit>`
 - :download:`Si_corehole.PBE-paw.corewf.abinit <_static/tutorial_14/Si_corehole.PBE-paw.corewf.abinit>`


References
----------

.. [1] N. D. M. Hine, Linear-Scaling Density Functional Theory using the Projector
  Augmented Wave Method, J. Phys. Condens. Matter 29, 024001 (2017).
  https://iopscience.iop.org/article/10.1088/0953-8984/29/2/024001
.. [2] Linear-Scaling Density Functional Theory and Theoretical Electron Energy
  Loss Spectroscopy Investigations of Surfaces and Defects in Nanomaterials,
  PhD Thesis of Edward Tait, 2019, University of Cambridge
  https://www.repository.cam.ac.uk/items/fcc71788-1468-47f7-989c-9a9d6e349f1e
