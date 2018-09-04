Example: SrO `IBRION 6` Calculation
===================================

This example uses the `OUTCARToGULP.py` script to process the `OUTCAR` file from an `IBRION 6 + ISIF 3` calculation on SrO.

Running the VASP calculation
----------------------------

The VASP `INCAR`, `KPOINTS` and `POSCAR` files used to perform the calculation are provided, and should be paired with the `Sr_sv` (`PAW_PBE Sr_sv 07Sep2000`) and `O` (`PAW_PBE O 08Apr2002`) PAW pseudopotentials from the VASP database.

The calculation is set up to perform a finite-difference phonon calculation (`IBRION = 6`) on a volume-relaxed optimised structure.
The displacement step `POTIM` is set to 0.01 &#8491; and `NFREE` is set to 4, so each atom will be displaced by &plusmn; 0.01 &#8491; and &plusmn; 0.02 &#8491; along the symmetry-unique directions.
`ISIF = 3` is also set, which instructs VASP to perform an additional set of lattice distortion to derive the elastic-constant matrix.

Running the calculation should produce a result similar to that in `OUTCAR.ref`.

Processing the output
---------------------

Running:

`python OUTCARToGULP.py -f OUTCAR.ref`

produces the following output:

```
Reading "OUTCAR.ref"...

Summary:
  -> Formula: SrO
  -> # structures: 33
  -> Stress tensors? Yes
  -> Phonon frequencies/eigenvectors? Yes
  -> Elastic-constant matrix? Yes

Writing data to "SrO.gulp"...

Done!
```

The data extracted from the `OUTCAR` file is written to `SrO.gulp`.

At the top of the file, a block of data is introduced by the comment:

```
# Calculated observables for SrO
# ==============================
```

This contains the "main" (initial) structure along with an `observables` block with the total energy, stress tensor, gradients, the unique non-zero components of the calculated elastic-constant matrix and the calculated phonon frequencies and eigenvectors.

The elastic constants are also used to calculate bulk and shear moduli using the Voigt definitions:

```
observables
    ...
    bulk_modulus
        94.24505
    shear_modulus
        60.64590
```

Below this are several blocks of data with headings similar to the following:

```
# Calculated energy, gradients and stresses for SrO (Structure 2)
# ===============================================================
```

These are the intermediate distorted/displaced structures along with the associated energies, gradients and stresses.

In this example, structures 2-25 correspond to lattice distortions, leading to changes in the stress tensor. Structures 26-33 correspond to displacement of the O and Sr atoms along the <i>x</i> direction, producing opposing forces (gradients) on both atoms.

Using optional command-line arguments to control the output
-----------------------------------------------------------

`OUTCARToGULP.py` has a number of additional options that can be used to customise the GULP input file, the usage of which is outlined below.

<h3>Adding basic commands</h3>

The `--add_commands` flag adds some basic GULP commands to the file:

`python OUTCARToGULP.py -f OUTCAR.ref -o SrO-Commands.gulp --add_commands`

With this option, the output `SrO-Commands.gulp` contains a `fit` command at the top, with the options based on the data it contains:

`fit conp comp prop phonon eigenvectors`

Some useful `output` commands are also added to the end of the file:

```
output thb fit.thb
output cif fit.cif
dump fit.grs
```

<h3>Adding fitting weights</h3>

For "hard" materials, the values of the elastic constants are likely to be large in comparison to other quantities such as stresses or gradients.

To avoid them being too heavily weighted during the fitting procedure, an optional weight can be added to the obserables block using the `--elastic_constants_weight` argument:

`python OUTCARToGULP.py -f OUTCAR.ref -o SrO-ECWeight.gulp --add_commands --elastic_constants_weight=1e-3`

This appends the supplied weight to the elastic constants in the observables block:

```
observables
    elastic 9
        1  1     186.09296
        1  2      48.32109
        ...
```

```
observables
    elastic 9
        1  1     186.09296  0.0010
        1  2      48.32109  0.0010
        ...
```

<h3>Limiting the output to the first structure</h3>

If the aim is to optimise a potential to reproduce the equilibrium structure and properties, the sequence of structures with lattice distortions and atomic displacements are probably of little interest.
These can be omitted using the `--first_structure` flag:

`python OUTCARToGULP.py -f OUTCAR.ref -o SrO-First.gulp --add_commands --first_structure`

The file `SrO-First.gulp` contains only one data block with the "main" structure and the calculated total energy, stress tensor, gradients, elastic-constant matrix and phonon frequencies/eigenvectors.

<h3>Limiting the output with gradient and stress thresholds</h3>

Alternatively, if the aim is to optimise a potential for lattice distortions or larger atomic displacements, it may be useful to extract only the configurations with large gradients or stresses.
To do this, the `--gradient_threshold` and `--stress_threshold` arguments can be used:

`python OUTCARToGULP.py -f OUTCAR.ref -o SrO-Threshold.gulp --add_commands --gradient_threshold=0.025 --stress_threshold=1.0`

The output file `SrO-Threshold.gulp` now contains data sets for a subset of the 33 structures, with the others replaced by single-line comments indicating they were omitted because the gradients/sstresses were below the set thresholds:

`# INFO: The gradient and/or stress-tensor components for structure 15 are below the set thresholds (gradients: 2.50e-02, stress: 1.00e+00) -> data set not output.`
