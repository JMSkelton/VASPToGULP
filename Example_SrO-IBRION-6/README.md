Example: SrO `IBRION 6` Calculation
===================================

This example uses the `OUTCARToGULP.py` script to process the `OUTCAR` file from an `IBRION 6` calculation on SrO.

Running the VASP calculation
----------------------------

The VASP `INCAR`, `KPOINTS` and `POSCAR` files used to perform the calculation are provided, and should be paired with the `Sr_sv` (`PAW_PBE Sr_sv 07Sep2000`) and `O` (`PAW_PBE O 08Apr2002`) PAW pseudopotentials from the VASP database.

The calculation is set up to perform a finite-difference phonon calculation (`IBRION = 6`) on a volume-relaxed optimised structure.
The displacement step `POTIM` is set to 0.01 &#8491; and `NFREE` is set to 4, so each atom will be displaced by &plusmn; 0.02 &#8491; and &plusmn; 0.01 &#8491; along the symmetry-unique directions.
`ISIF = 3` is also set, which requests VASP to perform an additional set of lattice distortions and to calculate the elastic-constant matrix.

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

This contains the "main" (initial) structure along with the unique non-zero components of the calculated elastic-constant matrix and the calculated phonon frequencies and eigenvectors in an `observables` block.

Below this are several blocks of data with headings similar to the following:

```
# Calculated gradients/strain derivatives for Structure 1
# =======================================================
```

These are the initial and intermediate displaced structures along with the associated gradients/strain derivatives, where these are above the default thresholds (10<sup>-5</sup> eV &#8491;<sup>-1</sup>).

In this example, structures 2-25 correspond to lattice distortions, leading to changes in the stress tensor (strain derivatives). Structures 26-33 correspond to displacement of the O and Sr atoms along the <i>x</i> direction, producing opposing forces (gradients) on both atoms.

Using the gradient/stress thresholds to control the output
----------------------------------------------------------

For fitting potentials, it is useful to extract only the configurations with large gradients/strain derivatives.

To do this, the `--gradient_threshold` and `--stress_threshold` command-line arguments can be used:

`python OUTCARToGULP.py -f OUTCAR.ref -o SrO-Threshold.gulp --gradient_threshold=0.025 --stress_threshold=0.1`

The output `SrO-Threshold.gulp` now contains data sets for a subset of the 33 structures, with the others replaced by single-line comments indicating they were omitted because the gradients/strain derivatives were below the set thresholds:

`# INFO: The gradient components and diagonal stress-tensor elements for structure 1 are below the set thresholds (gradients: 2.50e-02, stress: 1.00e-01) -> data set not output.`

Adding fitting weights
----------------------

For "hard" materials, the values of the elastic constants are likely to be large in comparison to other quantities such as stresses or gradients.

To avoid them being too heavily weighted during the fitting procedure, an optional weight can be added to the obserables block written to the output files using the `--elastic_constants_weight` command-line argument:

`python OUTCARToGULP.py -f OUTCAR.ref -o SrO-ECWeight.gulp --elastic_constants_weight=1e-3`

This appends weights to the elastic constants in the observables block:

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