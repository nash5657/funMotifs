#!/bin/bash
source /home/markmzr/.bashrc
source activate funMotif
echo "Start funMotifsMain"
module load PostgreSQL
python3 funMotifsMain.py --param_file='../conf/main_parameters_hg19.conf' --temp_dir='tmp' --force_overwrite=true >results/run14.txt 2>results/run14.err
