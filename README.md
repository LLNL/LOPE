# Lattice Optimization for Porous Electrodes (LOPE)

## Overview

This code performs an optimization over the structure of a porous electrode formed from a lattice of unit cells. This is a steepest descent version of the code used in:

V. A. Beck, J. J. Wong, C. F. Jekel, D. A. Tortorelli, S. E. Baker, E. B. Duoss and M. A. Worsley, Computational design of microarchitected porous electrodes for redox flow batteries, _J. Power Sources_, 2021, **512**, 230453. 


Please cite the reference above https://doi.org/10.1016/j.jpowsour.2021.230453, and/or the code https://doi.org/10.11578/dc.20220413.1 if using any of this code in your own work.

DOI: https://doi.org/10.11578/dc.20220413.1

## Build

The code requires installation of OpenFOAM. This code was tested using OpenFOAM v4.1. The code is compiled using the standard openFOAM command "wmake"

```
cd echemLatticeSmooth
wmake

```

## Run

The compiled code can be run after building a computational mesh

```
cd run01
blockMesh
echemLatticeSmoothSD

```

LLNL Release Number: LLNL-CODE-833809
