#!/bin/bash
#SBATCH -t 2:00:00
#SBATCH -N 1
#SBATCH -o /work2/10000/budvin/frontera/partition-project/paralab-partition-repo/mpi-cpp/output.txt
#SBATCH -e /work2/10000/budvin/frontera/partition-project/paralab-partition-repo/mpi-cpp/error.txt

#SBATCH -p development



#SBATCH --mail-user=u1472227@utah.edu  
#SBATCH --mail-type=FAIL




set -e

cd $SCRATCH/meshes
wget --no-check-certificate https://archive.nyu.edu/bitstream/2451/44304/1/10k_hex.tar.gz

tar -xvzf 10k_hex.tar.gz
