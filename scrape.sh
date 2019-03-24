#!/bin/env bash

#takes one argument -> projectname to search for (projname)-energy-run.out

#declaration of file inputs--list of simulation filenames, name of output file
filenames="filenames.txt"
simul_out="${1}-energy-run.out" #change on cluster
distances="distances.txt"
i=1

# each line is (trial num.)    (distance)     (energy)
while IFS= read -r line; do

  #write number and radial distance
  first=$(grep "^$i[^0-9]" $distances | tr -d '\n')
  echo -e -n "${first}\t" >> energies.txt

  #scrape energy from cp2k output and write
  cd $line
  grep "Total energy" $simul_out | grep -o -- "-[0-9]*\.[0-9]*" >> ../energies.txt

  #return to master folder
  cd ..
  i=$((i+1))

#for each line of submissions.txt
done < "$filenames"
