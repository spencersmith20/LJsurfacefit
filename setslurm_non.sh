#!/bin/env bash

slurm="cp2k-template.slurm"
inp="energy-geoopt.inp"
proj_name=${1}

for i in *.xyz; do

  #store the name (w/o .xyz) of the file in submissions.txt
  echo "${i%.*}" >> submissions.txt

  #create a directory of the .xyz file's name
  mkdir "${i%.*}"

  #move the .xyz file into it
  mv "${i}" "${i%.*}"

  #copy in the .inp file
  cp $slurm $inp "${i%.*}"

  #move to the created directory
  cd "${i%.*}"

  #replace inner arguments with sed
  sed -e "s|FIRST_PROJNAME|"${proj_name}"|g" \
      -e "s|EXTERNAL_XYZFILE|"${i}"|g" \
      -e "s|EXTERNAL_INPUTFILE|"energy-geoopt.inp"|g" \
      $slurm > "cp2k-submit.slurm"

  #submit the script to the cluster
  sbatch cp2k-submit.slurm

  #return to bulk directory
  cd ..

done
