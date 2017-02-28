# OUTCARToGULP.py by J. M. Skelton


# TODO: The _ReadOUTCARFile() function includes some basic sanity checks where it was convenient to do so, but these are certainly not exhaustive - it could probably include more validation.
# TODO: The _WriteGULPInputFile() function does not currently validate the dataSets parameter, which it ought to do.


import math;
import os;
import re;

from argparse import ArgumentParser;


# Regex to read the atom types from the summary of the POTCAR files written near the top of the OUTCAR file.

_VRHFINRegex = re.compile(r"VRHFIN =(?P<atom_type>[^\:]+)\:");

# Regex to extract the atomic masses from the summary of the POTCAR files.
# The semicolon in this expression is to avoid matching the line "POMASS = <mass_1>, <mass_2>, ..." that occurs later in the file.
# [In principle, however, the first parsing loop in _ReadOUTCARFile() should break before it hits this.]

_POMASSRegex = re.compile(r"POMASS =\s*(?P<pomass_value>\d+\.\d+);");

# Regex to read the NIONS tag printed near the top of the OUTCAR file.

_NIONSRegex = re.compile(r"NIONS =\s*(?P<nions_value>\d+)");

# Regex to read total energies.

_TOTENRegex = re.compile("TOTEN  =\s*(?P<toten_value>-?\d+\.\d+)\s+eV");

# Generic integer regex to capture the numbers of ions of each atom type, printed starting with "ions per type =".

_GenericIntegerRegex = re.compile(r"(?P<value>\d+)");

# Regex to match the lines in the OUTCAR file contaning the phonon frequencies.

_ModeFrequencyRegex = re.compile(r"^\s*(?P<mode_index>\d+) (?P<freq_sign>f|f/i)\s*=\s*(?P<freq_thz>\d+\.\d+) THz\s*(?P<freq_radthz>\d+\.\d+) 2PiTHz\s*(?P<freq_invcm>\d+\.\d+) cm-1\s*(?P<freq_mev>\d+\.\d+) meV\s*$");


# Function to parse a VASP OUTCAR file and extract information that can be imported into GULP.

def _ReadOUTCARFile(filePath = r"OUTCAR"):
    # Basic information.

    atomTypes = None;
    atomicMasses = None;
    nions, ionsPerType = None, None;

    # Stress tensors.

    stressTensors = None;

    # Structures (with sets of forces).

    structures = None;

    # Phonon frequencies and eigenvectors.

    frequencies, eigenvectors = None, None;

    # Elastic-constant matrix.

    elasticConstantMatrix = None;

    with open(filePath, 'r') as inputReader:
        # First, scan through the OUTCAR file and extract the atom types, atomic masses, the number of ions, and the number of ions of each type.
        # This information should be found close to the top of the file.

        for line in inputReader:
            if "VRHFIN" in line:
                match = _VRHFINRegex.search(line);

                if match:
                    if atomTypes == None:
                        atomTypes = [];

                    atomTypes.append(
                        match.group('atom_type').strip()
                        );
                else:
                    print("WARNING: _ReadOUTCARFile(): Failed to extract an atom type from a line containing the string 'VRHFIN'.");

            elif "POMASS" in line:
                match = _POMASSRegex.search(line);

                if match:
                    if atomicMasses == None:
                        atomicMasses = [];

                    atomicMasses.append(
                        float(match.group('pomass_value'))
                        );
                else:
                    print("WARNING: _ReadOUTCARFile(): Failed to extract a mass value from a line containing the string 'POMASS'.");

            elif "NIONS" in line:
                match = _NIONSRegex.search(line);

                if match:
                    if nions != None:
                        print("WATNING: _ReadOUTCARFile(): nions will be overwritten.")

                    nions = int(match.group('nions_value'));
                else:
                    print("WARNING: _ReadOUTCARFile(): Failed to extract a number of ions from a line containing the string 'NIONS'.");

            elif "ions per type =" in line:
                matches = _GenericIntegerRegex.findall(line);

                if len(matches) > 0:
                    if ionsPerType != None:
                        print("WARNING: _ReadOUTCARFile(): ionsPerType will be overwritten.");

                    ionsPerType = [int(match) for match in matches];

            # The number of atoms and the number of atoms per type should occur after the POTCAR summaries from which the atom types and atomic masses are extracted.
            # Once we have these, we can should have collected all the basic information we need, so we can move on.

            if atomTypes != None and atomicMasses != None and nions != None and ionsPerType != None:
                break;

        # Sanity checks.

        if sum(ionsPerType) != nions:
            raise Exception("Error: Inconsistent NIONS and \"ions per type\" lines in input file \"{0}\".".format(filePath));

        if len(atomTypes) != len(ionsPerType):
            raise Exception("Error: Inconsistent number of atom types read from input file \"{0}\".".format(filePath));

        if len(atomTypes) != len(ionsPerType):
            raise Exception("Error: Inconsistent number of atomic masses read from input file \"{0}\".".format(filePath));

        # Scan through the remainder of the file and capture useful blocks of information.

        for line in inputReader:
            line = line.strip();

            # The stress tensors are captured separately because certain values of the ISIF tag mean VASP does not calculate them, and they may then not appear in the output (I don't know - I've never tested this!).

            if line == "FORCE on cell =-STRESS in cart. coord.  units (eV):":
                # Block should contain a stress tensor.

                if stressTensors == None:
                    stressTensors = [];

                # Scan to the line which gives the total force, then store the stress tensor elements and break.

                for line in inputReader:
                    if line.strip()[:5] == "Total":
                        # VASP outputs the cell forces (= *minus* the stress), but the stress is more intuitive to work with -> we convert it here.

                        stressTensors.append(
                            [-1.0 * float(element) for element in line.strip().split()[1:7]]
                            );

                        break;

            if line == "VOLUME and BASIS-vectors are now :":
                # Block should contain a structure, a set of forces and a total energy.

                if structures == None:
                    structures = [];

                # Skip four lines.

                for i in range(0, 4):
                    next(inputReader);

                # The following three lines contain the real- and reciprocal-space lattice vectors.

                latticeVectors = [];

                for i in range(0, 3):
                    latticeVectors.append(
                        tuple(float(element) for element in next(inputReader).strip().split()[:3])
                        );

                # Scan forward until the header of the table that summarises the atomic positions and forces.

                for line in inputReader:
                    if "POSITION" in line:
                        break;

                # Skip one line.

                next(inputReader);

                # The following nions lines contain the positions of and the forces on each ion.

                atomPositions, forceSet = [], [];

                for i in range(0, nions):
                    elements = next(inputReader).strip().split();

                    atomPositions.append(
                        tuple(float(element) for element in elements[:3])
                        );

                    forceSet.append(
                        tuple(float(element) for element in elements[3:6])
                        );

                # Scan forward until the header of the block containing the total energy.

                for line in inputReader:
                    if line.strip() == "FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)":
                        break;

                # Skip one line.

                next(inputReader);

                # Extract the total energy.

                match = _TOTENRegex.search(next(inputReader));

                totalEnergy = float(match.group('toten_value'));

                # Store the data.

                structures.append(
                    (latticeVectors, atomPositions, totalEnergy, forceSet)
                    );

            elif line == "Eigenvectors and eigenvalues of the dynamical matrix":
                # Block contains frequencies and eigenvectors.

                # There should only be one set of frequencies and eigenvectors in a given OUTCAR file.

                if frequencies != None:
                    print("WARNING: _ReadOUTCARFile(): frequencies/eigenvectors will be overwritten.");

                frequencies, eigenvectors = [], [];

                # Skip three lines.

                for i in range(0, 3):
                    next(inputReader);

                # There should be 3 * nions eigenvectors to extract.

                for i in range(0, 3 * nions):
                    # The first line contains the phonon frequency in various units.

                    match = _ModeFrequencyRegex.match(next(inputReader));

                    if match:
                        # Sanity check.

                        if int(match.group('mode_index')) != i + 1:
                            raise Exception("Error: Unexpected mode index encountered while reading frequency/eigenvector block in \"{0}\".".format(filePath));

                        # We want the frequency in cm^-1.

                        frequency = float(match.group('freq_invcm'));

                        frequencies.append(
                            -1.0 * frequency if match.group('freq_sign') == 'f/i' else frequency
                            );
                    else:
                        raise Exception("Error: _ReadOUTCARFile(): Parsing phonon frequency line failed.");

                    # Skip one line.

                    next(inputReader);

                    # Read the eigenvector.

                    eigenvector = [];

                    # Each eigenvector should have nions components.

                    for j in range(0, nions):
                        eigenvector.append(
                            tuple(float(element) for element in next(inputReader).strip().split()[3:6])
                            );

                    eigenvectors.append(eigenvector);

                    # Skip the trailing blank line after the eigenvector.

                    next(inputReader);

            elif line == "TOTAL ELASTIC MODULI (kBar)":
                # Block contains the elastic-constant matrix.

                # Skip two lines.

                for i in range(0, 2):
                    next(inputReader);

                elasticConstantMatrix = [];

                for i in range(0, 6):
                    elasticConstantMatrix.append(
                        [float(element) for element in next(inputReader).strip().split()[1:7]]
                        );

        # If stress tensors were extracted, merge them with the structural data.
        # If not, add None values in their place.

        if stressTensors != None:
            if len(stressTensors) != len(structures):
                raise Exception("Error: Number of stress tensors does not match number of structures in input file \"{0}\".".format(filePath));

            for i, (latticeVectors, atomPositions, totalEnergy, forceSet) in enumerate(structures):
                structures[i] = (latticeVectors, atomPositions, totalEnergy, stressTensors[i], forceSet);
        else:
            for i, (latticeVectors, atomPositions, totalEnergy, forceSet) in enumerate(structures):
                structures[i] = (latticeVectors, atomPosition, totalEnergy, None, forceSet);

        # Generate a chemical formula for the structure.

        formula = "";

        for atomType, atomCount in zip(atomTypes, ionsPerType):
            if atomCount == 1:
                formula = formula + "{0}".format(atomType);
            else:
                formula = formula + "{0}{1}".format(atomType, atomCount);

        # Generate a list of atom types.

        atomTypesList = [];

        for atomType, atomCount in zip(atomTypes, ionsPerType):
            atomTypesList = atomTypesList + [atomType] * atomCount;

        # Generate a list of atomic masses.

        atomicMassesList = [];

        for atomicMass, atomCount in zip(atomicMasses, ionsPerType):
            atomicMassesList = atomicMassesList + [atomicMass] * atomCount;

        # Return the captured data.

        return (
            formula, atomTypesList, atomicMassesList,
            structures,
            (frequencies, eigenvectors) if frequencies != None else None,
            elasticConstantMatrix
            );

# Function to output sets of data in a GULP-friendly format.

# dataSets is a list of dictionaries containing some or all of the following keys:
#   -> 'HeaderComment' : a comment (required).
#   -> 'Name' : a name for the structure (required for output).
#   -> 'LatticeVectors' : 3 x three-component vectors (required for output).
#   -> 'AtomTypesList' : a list of atomic symbols (required for output).
#   -> 'AtomPositions' : a list of three-component vectors (required for output).
#   -> 'TotalEnergy' : the total energy.
#   -> 'StressTensor' : the unique elements of a stress tensor (xx, yy, zz, xy, yz, zx; optional).
#   -> 'ForceSet' : a list of three-component vectors (optional).
#   -> 'ElasticConstantMatrix' : a 6x6 matrix (optional).
#   -> 'PhononModes' : a tuple of (frequencies, eigenvectors) (optional).

# Each entry in the list must contain at least the 'HeaderComment' key.
# If 'Name', 'LatticeVectors', 'AtomTypes', and 'AtomPositions' are present, a minimal structure data block is output; if not, just the header comment is written.
# If any of the optional 'TotalEnergy', 'StressTensor', 'ForceSet', 'ElasticConstantMatrix' or 'PhononModes' keys are present, these are added to the data block as "observables".

# observablesWeights is an optional dictionary that can specify the weight assigned to a subset of the keys in dataSets.
# At present, only 'ElasticConstantMatrix' will be used, if present.

# addComments is a flag that adds some standard GULP comments to the output file when set.

# This code was refactored into a function to make it easier to arrange the output files in different ways depending on the data captured from the input file.

def _WriteGULPInputFile(dataSets, filePath, observablesWeights = None, addCommands = False):
    with open(filePath, 'w') as outputWriter:
        # If addComments is set, add a "fit" command to the top of the file.

        if addCommands:
            # The command options depend on what observables are to be written out.

            hasTotalEnergies, hasForceSets, hasStressTensorsOrElasticConstantMatrix, hasPhononModes = False, False, False, False;

            for dataSet in dataSets:
                if not hasTotalEnergies and 'TotalEnergy' in dataSet:
                    hasTotalEnergies = True;

                if not hasForceSets and 'ForceSet' in dataSet:
                    hasForceSets = True;

                if not hasStressTensorsOrElasticConstantMatrix and ('StressTensor' in dataSet or 'ElasticConstantMatrix' in dataSet):
                    hasStressTensorsOrElasticConstantMatrix = True;

                if not hasPhononModes and 'PhononModes' in dataSet:
                    hasPhononModes = True;

                if hasTotalEnergies and hasForceSets and hasStressTensorsOrElasticConstantMatrix and hasPhononModes:
                    break;

            # Build a list of fit options.

            fitKeywords = [];

            if hasForceSets or hasPhononModes:
                if hasStressTensorsOrElasticConstantMatrix:
                    fitKeywords.append('conp');
                else:
                    fitKeywords.append('conv');

                fitKeywords = fitKeywords + ['comp', 'prop'];

                if hasPhononModes:
                    fitKeywords = fitKeywords + ['phonon', 'eigenvectors'];

            if len(fitKeywords) == 0 and hasTotalEnergies:
                fitKeywords = ['fit', 'single', 'comp', 'prop'];

            # If there is something to fit, fitKeywords will not be empty -> write out the command.

            if len(fitKeywords) > 0:
                outputWriter.write(
                    "fit {0}\n".format(" ".join(fitKeywords))
                    );

                outputWriter.write("\n\n");
            else:
                # If not, reset addCommands and print a warning.

                print("WARNING: _WriteGULPInputFile(): addCommands was set but no observables were found in any of the items in dataSets.");

                addCommands = False;

        # To help format the output file, keep track of whether or not we are writing out the first data set, and whether the previous data set was a data block (or just a header comment in its place).

        firstDataSet, previousIsDataBlock = True, False;

        for dataSet in dataSets:
            isDataBlock = 'Name' in dataSet and 'LatticeVectors' in dataSet and 'AtomTypesList' in dataSet and 'AtomPositions' in dataSet;

            if not firstDataSet:
                if previousIsDataBlock or isDataBlock:
                    # We want data blocks to be visually separated from other data blocks/replacement comments, but data blocks replaced by comments can be bunched together.

                    outputWriter.write('\n\n');
            else:
                firstDataSet = False;

            # Output the header comment.

            outputWriter.write("# {0}\n".format(dataSet['HeaderComment']));

            if isDataBlock:
                # Underline the header comment and add a blank line below it.

                outputWriter.write("# {0}\n".format("=" * len(dataSet['HeaderComment'])));
                outputWriter.write('\n');

                # Output the data set name.

                outputWriter.write("name {0}\n".format(dataSet['Name']));
                outputWriter.write("\n");

                # Output the lattice vectors.

                outputWriter.write("vectors Angs\n");

                for x, y, z in dataSet['LatticeVectors']:
                    outputWriter.write("    {0: >15.9f}  {1: >15.9f}  {2: >15.9f}\n".format(x, y, z));

                outputWriter.write("\n");

                # Output the atom positions.

                outputWriter.write("cart\n");

                for atomType, (x, y, z) in zip(dataSet['AtomTypesList'], dataSet['AtomPositions']):
                    outputWriter.write("    {0: <3}  core  {1: >10.5f}  {2: >10.5f}  {3: >10.5f}\n".format(atomType, x, y, z));

                outputWriter.write("\n");

                # VASP doesn't print the actual spacegroup number, so we set this to 1 (P1) as a "safe" dummy value.

                outputWriter.write("space\n");
                outputWriter.write("1\n");

                # Check whether the data set contains "observables", and output if present.

                if 'StressTensor' in dataSet or 'ForceSet' in dataSet or 'ElasticConstantMatrix' in dataSet or 'PhononModes' in dataSet:
                    # Output a blank line to terminate the previous section.

                    outputWriter.write("\n");

                    outputWriter.write("observables\n");

                    # If available, output the total energy.

                    if 'TotalEnergy' in dataSet:
                        outputWriter.write("    energy eV\n");
                        outputWriter.write("        {0:.8f}\n".format(dataSet['TotalEnergy']));

                    # If available, output the diagonal components of the stress tensor.

                    if 'StressTensor' in dataSet:
                        outputWriter.write("    strain_derivative eV\n");

                        for j, element in enumerate(dataSet['StressTensor'][:3]):
                            outputWriter.write("        {0}  {1: >12.5f}\n".format(j + 1, element));

                    # If available, convert the forces to gradients and output.

                    if 'ForceSet' in dataSet:
                        outputWriter.write("    gradients eV/Angs\n");

                        for j, (fx, fy, fz) in enumerate(dataSet['ForceSet']):
                            # The gradients are *minus* the forces.

                            outputWriter.write("        {0: <3}  {1: >10.5f}  {2: >10.5f}  {3: >10.5f}\n".format(j + 1, -1.0 * fx, -1.0 * fy, -1.0 * fz));

                    # If available, output the non-zero elements of the upper triangle of the elastic-constant matrix.

                    if 'ElasticConstantMatrix' in dataSet:
                        ecMatrix = dataSet['ElasticConstantMatrix'];

                        data = [];

                        for i in range(0, 6):
                            for j in range(i, 6):
                                if ecMatrix[i][j] != 0.0:
                                    data.append(
                                        (i + 1, j + 1, ecMatrix[i][j])
                                        );

                        outputWriter.write("    elastic {0}\n".format(len(data)));

                        # If a weight for the elastic constants has been supplied via observablesWeights, add it.

                        ecWeight = None;

                        if observablesWeights != None and 'ElasticConstantMatrix' in observablesWeights:
                            ecWeight = observablesWeights['ElasticConstantMatrix'];

                        for i, j, cij in data:
                            # The elastic constants need to be in GPa rather than kbar -> divide by 10.

                            cij = cij / 10.0;

                            if ecWeight != None:
                                outputWriter.write("        {0}  {1}  {2: >12.5f}  {3:.4f}\n".format(i, j, cij, ecWeight));
                            else:
                                outputWriter.write("        {0}  {1}  {2: >12.5f}\n".format(i, j, cij));

                        # We also calculate the bulk and shear moduli and add them as observables.

                        c11, c22, c33, c44, c55, c66 = [ecMatrix[i][i] for i in range(0, 6)];
                        c12, c23, c31 = ecMatrix[0][1], ecMatrix[1][2], ecMatrix[2][0];

                        # The formulae used here are for the Voigt averages and were obtained from the Materials Project Wiki at https://materialsproject.org/wiki/index.php/Elasticity_calculations#Derived_elastic_properties.

                        bulkModulus = (1.0 / 9.0) * ((c11 + c22 + c33) + 2.0 * (c12 + c23 + c31)) / 10.0;

                        outputWriter.write("    bulk_modulus\n");
                        outputWriter.write("        {0:.5f}\n".format(bulkModulus));

                        shearModulus = (1.0 / 15.0) * ((c11 + c22 + c33) - (c12 + c23 + c31) + 3.0 * (c44 + c55 + c66)) / 10.0;

                        outputWriter.write("    shear_modulus\n");
                        outputWriter.write("        {0:.5f}\n".format(shearModulus));

                    # If available, output the phonon frequencies and eigenvectors.

                    if 'PhononModes' in dataSet:
                        frequencies, eigenvectors = dataSet['PhononModes'];

                        # Each phonon mode is input as a separate observable.

                        for frequency, eigenvector in zip(frequencies, eigenvectors):
                            outputWriter.write("    mode\n");

                            outputWriter.write("        {0: >12.6f}\n".format(frequency));

                            for dx, dy, dz in eigenvector:
                                outputWriter.write("            {0: >10.6f}  {1: >10.6f}  {2: >10.6f}\n".format(dx, dy, dz));

                    # Terminate the observables block.

                    outputWriter.write("end\n");

            # Update previousIsDataBlock.

            previousIsDataBlock = isDataBlock;

        # If addCommands is set, add some output commands to the end of the file.

        if addCommands:
            # Add two additional blank lines to terminate the previous section.

            outputWriter.write("\n\n");

            # Add output commands.

            outputWriter.write("output thb fit.thb\n");
            outputWriter.write("output cif fit.cif\n");
            outputWriter.write("dump fit.grs\n");


# Main.

if __name__ == "__main__":
    # Collect and parse command-line arguments.

    parser = ArgumentParser(description = "Extract fitting data for GULP from a VASP OUTCAR file and prepare a basic GULP input file");

    parser.set_defaults(
        InputFile = "OUTCAR",
        OutputFile = None,

        FirstStructure = False,
        AddCommands = False,

        GradientThreshold = None,
        StressThreshold = None,

        ElasticConstantsWeight = None
        );

    parser.add_argument(
        "-f", "--input_file",
        metavar = "<input_file>",
        type = str, dest = 'InputFile',
        help = "Input file to read (default: OUTCAR)"
        );

    parser.add_argument(
        "-o", "--output_file",
        metavar = "<output_file>",
        type = str, dest = 'OutputFile',
        help = "Output file (default: automatically determined)"
        );

    parser.add_argument(
        "--first_structure",
        action = 'store_true', dest = 'FirstStructure',
        help = "For OUTCAR files containing data for multiple structures, output data only for the first one (default: no)"
        );

    parser.add_argument(
        "--add_commands",
        action = 'store_true', dest = 'AddCommands',
        help = "Add some basic commands to the GULP input file (default: no)"
        );

    parser.add_argument(
        "--gradient_threshold",
        metavar = "<threshold>",
        type = float, dest = 'GradientThreshold',
        help = "Threshold for the output of forces (gradients); sets of gradients where all components have absolute values less than this will not be output (default: no threshold)"
        );

    parser.add_argument(
        "--stress_threshold",
        metavar = "<threshold>",
        type = float, dest = 'StressThreshold',
        help = "Threshold for the output of stress tensor components (strain derivatives); stress tensors where all components have absolute values less than this will not be output (default: no threshold)"
        );

    parser.add_argument(
        "--elastic_constants_weight",
        metavar = "<weight>",
        type = float, dest = 'ElasticConstantsWeight',
        help = "Fitting weight added to elastic constants observables in the output file (if written; default: no weight)"
        );

    args = parser.parse_args();

    # Perform some basic validation.

    if args.GradientThreshold != None and args.GradientThreshold < 0.0:
        raise Exception("Error: If supplied, the gradient component threshold must be >= 0.");

    if args.StressThreshold != None and args.StressThreshold < 0.0:
        raise Exception("Error: If supplied, the stress-tensor component threshold must be >= 0.");

    if args.ElasticConstantsWeight != None and args.ElasticConstantsWeight < 0.0:
        raise Exception("Error: If supplied, the elastic-constants weight must be >= 0.");

    # Read input file.

    print("Reading \"{0}\"...".format(args.InputFile));

    formula, atomTypesList, atomicMassesList, structures, phononModes, elasticConstantMatrix = _ReadOUTCARFile(args.InputFile);

    print("");

    # Print some information about what data was extracted from it.

    print("Summary:")
    print("  -> Formula: {0}".format(formula));
    print("  -> # structures: {0}".format(len(structures)));

    _, _, _, stressTensor, _ = structures[0];

    print("  -> Stress tensors? {0}".format("Yes" if stressTensor != None else "No"));

    print("  -> Phonon frequencies/eigenvectors? {0}".format("Yes" if phononModes != None else "No"));
    print("  -> Elastic-constant matrix? {0}".format("Yes" if elasticConstantMatrix != None else "No"));

    print("");

    # Build data sets to output.
    # The "rules" for doing this are loosely based on the various types of VASP calculation that are likely to be used to prepare data for fitting in GULP.

    dataSets = [];

    # If the --first_structure command-line argument was set, discard any additional structures read from the input file.

    if args.FirstStructure:
        structures = structures[:1];

    for i, (latticeVectors, atomPositions, totalEnergy, stressTensor, forceSet) in enumerate(structures):
        # Work out whether to output the gradients.

        outputGradients = False;

        # If no threshold is set, the gradients are always output.

        if args.GradientThreshold == None:
            outputGradients = True;
        else:
            # If a threshold has been set, the gradients are output if the absolute value of any component of any of the set is larger than the threshold.

            for fx, fy, fz in forceSet:
                if math.fabs(fx) >= args.GradientThreshold or math.fabs(fy) >= args.GradientThreshold or math.fabs(fz) >= args.GradientThreshold:
                    outputGradients = True;
                    break;

        # Work out whether to output the stress tensor.

        outputStressTensor = False;

        if stressTensor != None:
            # As for gradients, if no threshold is set, the diagonal components of the stress tensor are always output.

            if args.StressThreshold == None:
                outputStressTensor = True;
            else:
                # If a threshold has been set, the diagonal components of the stress tensor are output when any are above the threshold.

                sXX, sYY, sZZ, _, _, _ = stressTensor;
                outputStressTensor = math.fabs(sXX) >= args.StressThreshold or math.fabs(sYY) >= args.StressThreshold or math.fabs(sZZ) >= args.StressThreshold;

        # Determine whether to output the structure.

        # If we are outputting the gradients or stress tensor, output the structure.

        outputStructure = outputGradients or outputStressTensor;

        # If we have elastic constants and/or phonon frequencies/eigenvectors, we assume these properties were calculated for the first structure -> make sure it gets output.

        outputStructure = outputStructure or (i == 0 and (phononModes != None or elasticConstantMatrix != None));

        # If this is the only structure read from the file, make sure it gets output.

        outputStructure = outputStructure or len(structures) == 1;

        if outputStructure:
            # Devise a header comment and a name.

            name, headerComment = None, None;

            if i == 0 and (phononModes != None or elasticConstantMatrix != None):
                name = formula;
                headerComment = "Calculated observables for {0}".format(formula);
            elif i == 0 and len(structures) == 1:
                name = formula;
                headerComment = formula;
            else:
                name = "{0} (Structure {1})".format(formula, i + 1);
                headerComment = "Calculated energy, gradients and strain derivatives for {0}".format(name);

            dataSet = {
                'HeaderComment' : headerComment,
                'Name' : name,

                'LatticeVectors' : latticeVectors,
                'AtomTypesList' : atomTypesList,
                'AtomPositions' : atomPositions,

                'TotalEnergy' : totalEnergy
                };

            if outputStressTensor:
                dataSet['StressTensor'] = stressTensor;

            if outputGradients:
                dataSet['ForceSet'] = forceSet;

            if i == 0:
                if elasticConstantMatrix != None:
                    dataSet['ElasticConstantMatrix'] = elasticConstantMatrix;

                if phononModes != None:
                    dataSet['PhononModes'] = phononModes;

            dataSets.append(dataSet);
        else:
            # If the gradients and diagonal stress-tensor elements are below the threshold, output a comment to note why the data set was excluded.

            dataSets.append(
                { 'HeaderComment' : "INFO: The gradient components and diagonal stress-tensor elements for structure {0} are below the set thresholds (gradients: {1:.2e}, stress: {2:.2e}) -> data set not output.".format(i + 1, args.GradientThreshold, args.StressThreshold) }
                );

    # Work out a name for the output file.

    # If an output file has been specified through the command line, use that.

    outputFile = args.OutputFile;

    # If not, determine one based on the chemical formula extracted from the OUTCAR file.

    if outputFile == None:
        outputFile = "{0}.gulp".format(formula);

        # Append numbers to the file name to overwriting other output files.

        if os.path.isfile(outputFile):
            fileNumber = 2;

            while True:
                outputFile = "{0}-{1}.gulp".format(formula, fileNumber);

                if not os.path.isfile(outputFile):
                    break;

                fileNumber = fileNumber + 1;

    # Compile a set of observables weights, depending on the supplied command-line options.

    observablesWeights = { };

    if args.ElasticConstantsWeight != None:
        observablesWeights['ElasticConstantMatrix'] = args.ElasticConstantsWeight;

    # Write the output file.

    print("Writing data to \"{0}\"...".format(outputFile));

    _WriteGULPInputFile(dataSets, outputFile, observablesWeights = observablesWeights, addCommands = args.AddCommands);

    print("");

    # Print a "finished" message.

    print("Done!");
