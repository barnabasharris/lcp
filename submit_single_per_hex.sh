#!/bin/bash -l

#  R MPI parallel job

# Request wallclock time (format hours:minutes:seconds).
#$ -l h_rt=3:00:0

# Request 1 gigabyte of RAM per process.
#$ -l mem=2000M

# Set the name of the job.
#$ -N lcp

# Select the MPI parallel environment with 32 processes
#$ -pe mpi 35

# Load R/GRASS environment
echo "running init.sh script..."
source /home/tcrnbgh/Scratch/lcp/init.sh
echo "done!"

# Set the working directory to somewhere in your scratch space.  This is
# necessary because the compute nodes cannot write to your $HOME
# NOTE: this directory must exist.
# set working dir
#$ -wd /home/tcrnbgh/Scratch/lcp

# Run our MPI job. GERun is our wrapper for mpirun, which launches MPI jobs  
echo "running gerun..."
gerun /home/tcrnbgh/RMPISNOW_bgh < /home/tcrnbgh/Scratch/lcp/rscript/analysis_lcp_single_thread_per_hex_rmpi.R
echo "done!"