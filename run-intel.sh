#!/bin/bash

## Project:
#SBATCH --account=nn9516k

## QOS?
# #SBATCH --qos=devel

## Job name:
#SBATCH --job-name=LeeroyJenkins

## Wall time limit:
#SBATCH --time=0-23:00:00

## Number of nodes:
#SBATCH --nodes=10

## Number of tasks to start on each node:
#SBATCH --ntasks-per-node=1

## Set OMP_NUM_THREADS
#SBATCH --cpus-per-task=32

# #Memory per core
# #SBATCH --mem-per-cpu=199MB

# Mail notifications
#SBATCH --mail-type=ALL

## Recommended safety settings:
set -o errexit # Make bash exit on any error
set -o nounset # Treat unset variables as errors

## Software modules
module restore system   # Restore loaded modules to the default
ml load intel/2018b
module list             # List loaded modules, for easier debugging

export OMP_PROC_BIND=true
export OMP_PLACES=cores

##
mkdir ${SLURM_JOBID}
mkdir -p ${SCRATCH}
cd ${SCRATCH}

## Copy source files
cp ${SUBMITDIR}/pbucketmpi.f90 ${SCRATCH}
cp ${SUBMITDIR}/functionsmpi.f90 ${SCRATCH}
cp ${SUBMITDIR}/chaosbucketmpi.f90 ${SCRATCH}

## Compile code
mpiifort -no-wrap-margin -qopenmp -O3 pbucketmpi.f90 functionsmpi.f90 chaosbucketmpi.f90 -o bucket.out -lmkl_rt
# mpiifort -no-wrap-margin -qopenmp -g -check all -fpe0 -warn -traceback -debug extended pbucketmpi.f90 functionsmpi.f90 chaosbucketmpi.f90 -o bucket.out -lmkl_rt

## Run the application
# srun -n 1 bucket.out
mpirun ./bucket.out

## Copy data and source back to submitdir
cp *.f90 ${SUBMITDIR}/${SLURM_JOBID}
cp *.txt ${SUBMITDIR}/${SLURM_JOBID}

## Tidy up, backup source code
# tidy?

## Happy end
exit 0
