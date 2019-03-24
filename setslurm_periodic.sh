#!/bin/env bash

#input files
slurm="cp2k-per-template.slurm"
inp="periodic.inp"
box_inp="box.txt"
proj_name=${1}
box_vals=( `cat $box_inp `)

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
  sed -e "s|FIRST_PROJNAME|${proj_name}|g" \
      -e "s/EXTERNAL_XYZFILE/$i/g" \
      -e "s/EXTERNAL_INPUTFILE/$inp/g" \
      -e "s/EXTERNAL_BBOX_X/${box_vals[0]}/g" \
      -e "s/EXTERNAL_BBOX_Y/${box_vals[1]}/g" \
      -e "s/EXTERNAL_BBOX_Z/${box_vals[2]}/g" \
      $slurm > "cp2k-per-submit.slurm"

  #submit the script to the cluster
  sbatch cp2k-per-submit.slurm

  #return to bulk directory
  cd ..

done
