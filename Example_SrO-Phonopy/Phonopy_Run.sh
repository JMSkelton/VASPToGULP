#!/bin/bash

for step in 0.01 0.05 0.10 0.50 1.00
do
    phonopyDir="Step-${step}"
    
    if [ -d "${phonopyDir}" ] ; then
        cd "${phonopyDir}"
        
        for poscarFile in `ls POSCAR-*`
        do
            runDir="${poscarFile/POSCAR-/}"
            
            if [ ! -d "${runDir}" ] ; then
                mkdir "${runDir}"
                
                cd "${runDir}"
                
                cp "../${poscarFile}" "POSCAR"
                
                for inputFile in "INCAR" "KPOINTS" "POTCAR"
                do
                    cp "../../${inputFile}" .
                done
                
                mpirun -np 16 vasp_std | tee "../${runDir}.out"
                
                cd ..
            fi
        done
        
        cd ..
    fi
done

