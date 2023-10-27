#!/bin/bash
source /home/markmzr/.bashrc
source activate funMotif
module load PostgreSQL
pg_ctl restart -D /proj/snic2020-16-187/private/database
echo "Start funMotifsMain"
python3 funMotifsMain.py -m --param_file='../conf/main_parameters_MM2.conf' --temp_dir='tmp' >results/MA1.txt 2>results/MA1.err
pg_ctl stop -D /proj/snic2020-16-187/private/database
