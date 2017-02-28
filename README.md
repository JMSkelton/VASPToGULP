VASPToGULP
==========

A collection of scripts for converting output from the Vienna *ab initio* Simulation Package (VASP) to input for the General Utility Lattice Program (GULP).

`OUTCARToGULP.py`
-----------------

Reads a VASP `OUTCAR` file, extracts various sets of data, and writes them to a GULP-friendly output format.

This script currently collects the following data:

- Structures with lattice vectors, atom types/positions, total energies, forces and stress tensors.
- Elastic-constant matrices.
- Phonon modes (frequencies and eigenvectors).

For single-point calculations, the script writes out the structure plus the total energy, gradients (minus forces) and strain derivatives (diagonal elements of the stress tensor).

For finite-differences calculations (e.g. `IBRION = 6` with `ISIF = 3`), the script writes out the initial structure with the calculated elastic-constant matrix and/or phonon frequencies/eigenvectors, as well as each each of the intermediate structures with lattice distortions and atomic displacements.

Several optional command-line arguments allow the output to be customised - type `python OUTCARToGULP.py -h` for details, or look at the SrO/`IBRION 6` [example](./Example_SrO-IBRION-6), which illustrates the basic usage of the script.

`OUTCARToGULP_ModeMap.py`
-------------------------

An auxiliary script for processing the VASP `OUTCAR` files generated from a sequence of single-point calculations for mapping a phonon-mode potential-energy surface using the [ModeMap](https://github.com/JMSkelton/ModeMap) tool.

The script extracts the total energies, forces (gradients) and stress tensors (strain derivatives) from each structure and collects them into a single GULP-friendly output file.

It also accepts some of the same command-line arguments as `OUTCARToGULP.py`, which can be used to customise the output.

The usage of this script is illustrated by the [ModeMap example](./Example_SrO-ModeMap).
