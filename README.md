# JDAM
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4049976.svg)](https://doi.org/10.5281/zenodo.4049976)

Copyright (c) 2020, Xy Wang and Nadav Avidor.
All rights reserved.
JDAM - Jump Diffusion by Analytical Models, subject to the GNU/GPL-3.0-or-later.

Please cite JDAM software using the appropriate DOI.

JDAM is a Matlab based software that computes the Intermediate Scattering Function(ISF) of adsorbate surface jump diffusion processes.
The software calculates the lineshape for arbitrary lattices which can be defined by the user. At the moment, the package includes an example definition for a hexagonal surface, including the principle adsorption sites; top, hollow (hcp+fcc), and three bridge sites. Furthermore, at the moment JDAM computes jumps on that surface only to nearest neighbours. In addition, an explicit case for jumps on hollow sites only, up to 10th nearest neighbours, is provided.


####################
# Getting started: #
####################

Installation:
-------------
- Install Matlab (tested with 2019b).

Usage:
------------------
- Prepare a file to describe the lattice. An example is given in hex_ui.m
- Edit the user interface file, <lattice>_ui.m .
- Edit the run file, run_jdam.m, to determine both a few parameters, and the user interface file, as well as the file which describe the lattice.
- Run run_jdam.m

Visualization:
--------------
The program plots the lineshape as well as the decay rates and intensities of various components.

####################
# List of Files:   #
####################

hex_ui.m
hex_vectors.m
run_jdam.m
surf_gen.m
calc_ISF.m
multiple_hollow_A.m
plot_isf.m
