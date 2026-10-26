#!/bin/bash
#SBATCH --job-name=composite_t2
#SBATCH --partition=C064M0256G
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --output=/lustre/home/2101110354/sequence_entropy/project/mouse_preimp/mutual_information/logs/composite_t2_%j.out
#SBATCH --error=/lustre/home/2101110354/sequence_entropy/project/mouse_preimp/mutual_information/logs/composite_t2_%j.err
/lustre/home/2101110354/anaconda3/envs/hoi_env/bin/python3 /lustre/home/2101110354/sequence_entropy/project/mouse_preimp/mutual_information/scripts/plot_composite_trend_hist_n3_n5_n10_transposed.py
