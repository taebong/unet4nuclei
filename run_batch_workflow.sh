#!/bin/bash
#SBATCH -c 1
#SBATCH -N 1
#SBATCH -t 0-12:00:00
#SBATCH -p gpu
#SBATCH --gres=gpu:1
#SBATCH --mem=4000
#SBATCH -o log/batch_workflow_%j.out
#SBTACH -e log/batch_workflow_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=taeyeon_yoo@hms.harvard.edu

module purge
module load gcc/6.2.0
module load cuda/9.0
module load python/3.6.0

source nucleiseg_env/bin/activate

source batch_workflow.sh
