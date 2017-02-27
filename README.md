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

It is also possible to have a fitting weight added to the elastic constants using the `--elastic_constants_weight` argument.

The SrO/`IBRION 6` [example](./Example_SrO-IBRION-6) illustrates the basic usage of the script.

`OUTCARToGULP_ModeMap.py`
-------------------------

An auxiliary script for processing the VASP `OUTCAR` files generated from a sequence of single-point calculations for mapping a phonon-mode potential-energy surface using the [ModeMap](https://github.com/JMSkelton/ModeMap) tool.

The script extracts the forces (gradients) and stress tensors (strain derivatives) from each structure and collects them into a single GULP-friendly output file.

It also accepts the same `--gradient_threshold` and `--stress_threshold` command-line arguments as `OUTCARToGULP.py`, which can be used to select structures with large gradients and/or stresses if desired.

The usage of this script is best illustrated by the [ModeMap example](./Example_SrO-ModeMap).
