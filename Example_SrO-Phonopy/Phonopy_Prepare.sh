#!/bin/bash

for step in 0.01 0.05 0.10 0.50 1.00
do
    dir="Step-${step}"
    
    mkdir "${dir}"
    
    cd "${dir}"
    
    phonopy -d --dim="2 2 2" -c ../POSCAR.vasp --amplitude=${step}
    
    cd ..
done
