
========================================================
Tutorial 14: Density Functional Tight Binding in ONETEP
========================================================

:Author: Arihant Bhandari
:Date:   July 2024

.. role:: raw-latex(raw)
   :format: latex

Introduction
============

Density Functional Tight binding models have been derived by a Taylor expansion
of the DFT energy functional in terms of the electron density and truncation
up to a certain order in the expansion [Foulkes1989]_. Within DFTB, 
the following eigenvalue equations are solved via diagonalization:

.. math::
  :label: Kohn_Sham_equation

   \begin{aligned}
      H_{\alpha\beta}M^{\beta}_{\,i} = S_{\alpha\beta} M^{\beta}_{\,i} \epsilon_i,  
   \end{aligned}

where :math:`H_{\alpha\beta}` is the hamiltonian matrix, :math:`S_{\alpha\beta}` 
is the overlap matrix, :math:`M^{\beta}_{\,i}` are the orbital coefficients, and 
:math:`\epsilon_i` are energy eigenvalues. The Hamiltonian is built from
parameters. [Elstner1998]_ proposed a self consistent charge (SCC) extension 
to the traditional DFTB approach, which optmizes the atomic charges
self consistently. Henceforth, a series of SCC DFTB methods have been developed [Gaus2014]_. 
Recently, non-SCC methods have undergone a revival because of their speed and applicability 
to large systems [Bannwarth2020]_. E.g. GFN0 is one such method  
where atomic charges are found using a charge equilibriation scheme [Pracht2019]_. 
The total energy in GFN0 also includes zeroth order terms such as dispersion, repulsion, 
electrostatic interactions and short range basis correction. 

We have implemented the GFN0 method within ONETEP. 
Here we include D2 dispersion correction [Grimme2006]_ instead of D4 [Caldeweyher2019]_.
The standard ensemble-DFT subroutines are used for diagonalization and
calculation of electronic energies and forces. 

Keywords
========

Here are some basic keywords to perform a DFTB calculation.

-  ``dftb: T/F`` 
  | [Boolean, default ``dftb: F``]. 
  | If true, it enables DFTB calculations.

-  ``dftb_method: GFN0`` 
  | [Text, default ``dftb_method: GFN0``]. 
  | Variant of the DFTB method, only GFN0 has been implemented at the moment. 

-  ``dftb_method_param_file: file address`` 
  | [Text, default ``dftb_method_param_file: param_gfn0-xtb.txt``]. 
  | Path to the parameter file. A specimen file is supplied in utils-devel/dftb folder. 

-  ``dftb_common_param_file: file address`` 
  | [Text, default ``dftb_common_param_file: param_gfn_common.txt``]. 
  | Path to the file for common GFN parameters. A specimen file is supplied in utils-devel/dftb folder. 

-  ``dftb_bc: O O O / P P P`` 
  | [Boolean, default ``dftb_bc: P P P``]. 
  | Boundary conditions. Only fully open (O O O) or full periodic (P P P) are implemented. 


Input files
===========

In this tutorial, we will use DFTB to calculate the relaxed geometry of ethylene carbonate molecule (OBC),
lithium tetra ethylene carbonate cluster (OBC), and perform molecular dynamics on a lithium graphite system (PBC).
Please download the files below and run them using ONETEP. 

- :download:`ethylene_carbonate.dat <_static/tutorial_DFTB/ethylene_carbonate.dat>`
- :download:`liec4.dat <_static/tutorial_DFTB/liec4.dat>`
- :download:`li-graphite.dat <_static/tutorial_DFTB/li-graphite.dat>`

The DFTB-GFN0 parameter files are available with [utils-devel]_ repository and also below: 

- :download:`param_gfn0-xtb.txt <_static/tutorial_DFTB/param_gfn0-xtb.txt>`
- :download:`param_gfn_common.txt <_static/tutorial_DFTB/param_gfn_common.txt>`

References
==========

.. [Foulkes1989] \ W. Matthew C. Foulkes, Roger Haydock, *Phys. Rev. B* **1989**, 39, 12520, https://doi.org/10.1103/PhysRevB.39.12520

.. [Elstner1998] Marcus Elstner et. al., *Phys. Rev. B* **1998**, 58, 7260, https://doi.org/10.1103/PhysRevB.58.7260

.. [Gaus2014] Michael Gaus, Qiang Cui, Marcus Elstner, "Density functional tight binding: application to organic biological molecules", *WIREs Comput. Mol. Sci.* **2014**, 4, 49, https://doi.org/10.1002/wcms.1156

.. [Bannwarth2020] Christoph Bannwarth et. al., "Extended tight-binding quantum chemistry methods", *WIREs Comput. Mol. Sci.* **2021**, 11, 1, https://doi.org/10.1002/wcms.1493

.. [Pracht2019] Philipp Pracht, Eike Caldeweyher, Sebastian Ehlert, Stefan Grimme, "A robust non-self-consistent tight-binding quantum chemistry method for large molecules", *ChemRxiv* **2019**, https://doi.org/10.26434/chemrxiv.8326202.v1

.. [Grimme2006] Stefan Grimme, "Semi-empirical GGA-type density functional constructed with a long-range dispersion correction", *J. Comput. Chem.* **2006**, 27, 1787, https://doi.org/10.1002/jcc.20495

.. [Caldeweyher2019] Eike Caldeweyher et. al., "A generally applicable atomic-charge dependent London dispersion correction", *J. Chem. Phys.* **2019**, 150, 154122, https://doi.org/10.1063/1.5090222

.. [utils-devel] https://github.com/onetep-devel/utils-devel 

