#!/bin/env bash
#SBATCH --job-name=FIRST_PROJNAME-energy
#SBATCH --output=FIRST_PROJNAME-energy.out
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=18
#SBATCH --cluster=mpi
#SBATCH --time=0-04:00:00

## ======================================================
## ===== Only make changes here!!!!  ====================
## ======================================================
   ! Can also change job and output names above
runtype=ENERGY
sysname=FIRST_PROJNAME
xyzfile=EXTERNAL_XYZFILE
template=EXTERNAL_INPUTFILE
box_x=EXTERNAL_BBOX_X
box_y=EXTERNAL_BBOX_Y
box_z=EXTERNAL_BBOX_Z

project_file="${sysname}-energy"
print_level="MEDIUM"
   ! rename project file as needed
   ! use print_level "LOW" for geoopt and print_level "MEDIUM" or energy
charge="0"
lsd="FALSE"

## ======================================================
## ===== Prepare the Cores ==============================
## ======================================================
module purge
module load intel/2017.3.196 intel-mpi/2017.3.196 cp2k/6.1
ulimit -s unlimited
ulimit -l unlimited

## ======================================================
## ===== Copy all necessary files to all nodes ==========
## ======================================================
files=(${template} ${xyzfile})
for i in ${files[@]}; do
    sbcast $SLURM_SUBMIT_DIR/$i $SLURM_SCRATCH/$i
done

## ======================================================
## ===== Copy all files back at end of run ==============
## ======================================================
run_on_exit(){
    cp -r $SLURM_SCRATCH/* $SLURM_SUBMIT_DIR
    pkill --uid=$SLURM_JOB_USER cp2k.popt
}
trap run_on_exit EXIT

## ======================================================
## ===== Geometry Optimization Calculations =============
## ======================================================

cd $SLURM_SCRATCH
export I_MPI_FABRICS_LIST=ofa,tmi
export I_MPI_FALLBACK=0

## Fill in template
sed -e "s|EXTERNAL_XYZ_FILE|${xyzfile}|g" \
    -e "s/EXTERNAL_PROJNAME/${project_file}/g" \
    -e "s/EXTERNAL_RUN_TYPE/${runtype}/g" \
    -e "s/EXTERNAL_PRINT_LEVEL/${print_level}/g" \
    -e "s/EXTERNAL_CHARGE/${charge}/g" \
    -e "s/EXTERNAL_LSD/${lsd}/g" \
    -e "s/EXTERNAL_BOX_X/${box_x}|g" \
    -e "s/EXTERNAL_BOX_Y/${box_y}|g" \
    -e "s/EXTERNAL_BOX_Z/${box_z}|g" \
    ${template} > "${project_file}-run.inp"

## Run calculation
echo "Beginning ${runtype} Calculation"
mpirun -np $SLURM_NTASKS cp2k.popt -i ${project_file}-run.inp -o ${project_file}-run.out
echo "Calculation finished successfully"
