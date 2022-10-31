#!/bin/bash
source /home/markmzr/.bashrc
source activate funMotif
echo "Start funMotifsMain"
module load PostgreSQL
python3 funMotifsMain.py --param_file='../conf/main_parameters.conf' --temp_dir='tmp' --force_overwrite=true >results/run24.txt 2>results/run24.err
