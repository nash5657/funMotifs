#!/bin/bash
source /home/markmzr/.bashrc
source activate funMotif
module load PostgreSQL
pg_ctl restart -D /proj/snic2020-16-187/private/database
echo "Start funMotifsMain"
python3 funMotifsMain.py -r -f --param_file='../conf/main_parameters_MM.conf' --temp_dir='tmp' >results/test2.txt 2>results/test2.err
pg_ctl stop -D /proj/snic2020-16-187/private/database
