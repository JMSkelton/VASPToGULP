#!/bin/bash

for step in 0.01 0.05 0.10 0.50 1.00
do
    dir="Step-${step}"

    cd "${dir}"

    phonopy -f 001/vasprun.xml 002/vasprun.xml

    phonopy -c ../POSCAR.vasp --dim="2 2 2" --mesh="2 2 2" --gamma_center --fc_symmetry=1
    mv mesh.yaml mesh-2x2x2.yaml

    phonopy -c ../POSCAR.vasp -p -s ../Phonopy_PhononDoS.conf
    phonopy -c ../POSCAR.vasp -p -s ../Phonopy_PhononDispersion.conf

    cd ..
done
