FILE:					USE:
cp2k-per-template.slurm			modified template for periodic SPE calcs
cp2k-template.slurm			modified template for non-periodic SPE calcs
energy-geoopt.inp			non-periodic CP2K parameters file
filecreate_non.py			generates n concatenated xyz files periodic) with randomized molecular orientation in a specified window of z values, .txt on box dimensions, .txt on each 	numbered system name
filecreate_periodic.py			generates n concatenated xyz files (non-periodic) with randomized molecular orientation in a specified window of z values, .txt on radial distances, .txt on 	each numbered system name
fit.py					temporary least squares solver for energies
graphene.cif				cif file used for periodic calcs
graphene.xyz				graphene input file used for non-periodic testing
leastsquare.py				fitting script that calculates energies between the surface and the m bodies in the molecule. refines parameters using a nonlinear least sqaures algorithm 
periodic.inp				periodic CP2K parameters file
scrape.sh				scrapes energies from each energy .out file to create energies.txt file with trial number, radial distance, and SPE result **REQUIRES "PROJNAME" (see below)
setslurm_non.sh				move each generated .xyz file into system folder with needed input files and submit job to cluster **REQUIRES DESIRED "PROJNAME" i.e. "sotolon" for a cp2k result of "sotolon-energy-run.out"
setslurm_periodic.sh			same as setslurm_non.sh for periodic system (different filenames)
sotolon.xyz				sotolon xyz file used for testing

PROCEDURE:
-use one of the filecreate scripts to generate required .xyz files
-change into newly-created directory to use appropriate setslurm script to submit calculations i.e. bash setslurm_non.sh "sotolon"
-once jobs are completed, scrape.sh can be run from project folder to generate energies.txt
-use fitting script to generate LJ parameters for m-body molecule

updated 3/24/2019 by spencer smith
