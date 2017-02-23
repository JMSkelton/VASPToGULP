VASPToGULP
==========

A collection of scripts for converting output from the Vienna *ab initio* Simulation Package (VASP) to input for the General Utility Lattice Program (GULP).

`OUTCARToGULP.py`
-----------------

Reads a VASP `OUTCAR` file, extracts various sets of data, and writes them to a GULP-friendly output format.

This script currently collects the following data:

- Structures with lattice vectors, atom types/positions, forces and stress tensors.
- Elastic-constant matrices.
- Phonon modes (frequencies and eigenvectors).

For single-point calculations, the script writes out the structure, plus the gradients (minus forces) and/or strain derivatives (diagonal elements of the stress tensor) if these are non-zero.

For finite-differences calculations (e.g. `IBRION = 6` with `ISIF = 3`), the script writes out the initial structure with the calculated elastic-constant matrix and phonon frequencies/eigenvectors, followed by each structure for which the gradients and/or strain derivatives are non-zero.

Optional thresholds controlling when gradients/strain derivatives are output can be set via the `--gradient_threshold` and `--stress_threshold` command-line parameters.

The [SrO example](./Example_SrO-IBRION-6) illustrates the basic usage of the script.
