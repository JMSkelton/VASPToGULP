#!/bin/bash

for step in 0.01 0.05 0.10 0.50 1.00
do
    dir="Step-${step}"

    cd "${dir}"

    python OUTCARToGULP_ModeMap.py */OUTCAR -o "Step-${step}.gulp" -n "Phonopy (${step} A)" --add_commands

    cd ..
done

head -n 3 Step-0.01/Step-0.01.gulp > Phonopy.gulp

for step in 0.01 0.05 0.10 0.50 1.00
do
    tail -n +4 Step-${step}/Step-${step}.gulp | tail -r  | tail -n +4 | tail -r >> Phonopy.gulp

    # The version of `head` macOS 10.12.3 doesn't accept negative -n arguments; on Unix machines using versions that do, the following works and avoids the two calls to `tail -r`.
    # tail -n +4 Step-${step}/Step-${step}.gulp | head -n -3 >> Phonopy.gulp
done

tail -n 3 Step-0.01/Step-0.01.gulp >> Phonopy.gulp
