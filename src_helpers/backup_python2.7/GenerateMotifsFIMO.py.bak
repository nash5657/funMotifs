'''
Created on Feb 16, 2017

@author: husensofteng
'''
import sys
import os.path, subprocess
from subprocess import STDOUT,PIPE
import shlex
import numpy as np
from multiprocessing import Pool
import argparse

def get_sig_motif_instances(instances_file_fimo_format,
                            bed_input_file, 
                            sig_bed_output_file, 
                            sig_ranked_bed_output_file, 
                            limit_to_check, scores_sd_above_mean, 
                            percentage_highest_scored_isntances):
    '''Given a bed file of motif instances, this function returns a file containing 
    those that have a P-value significant, positive score and a score larger than: mean_scores+(sd_scores*scores_sd_above_mean.
    It also generates a file containing top x highest scored-instances based on the value of percentage_highest_scored_isntances'''
    score_index=6
    all_lines = []
    all_scores = []
    with open(instances_file_fimo_format, 'r') as bed_infile:
        l = bed_infile.readline()[1:]#skip the header line
        while l!= "":
            sl = l.strip().split('\t')
            if len(sl)==8:
                if float(sl[score_index])>0.0:#only include motifs that have a score above zero (ignore negative scores)
                    lw = []
                    for i in range(0, len(sl)):
                        if i == score_index:
                            lw.append(float(sl[i]))#append column i to list lw
                            all_scores.append(float(sl[i]))#if it was the 6th col (score) then append it to all_scoers 
                        else:
                            lw.append(sl[i])#append column i to list lw
                    all_lines.append(lw)
            l = bed_infile.readline()
    mean_scores = 0.0
    sd_scores = 0.0
    
    lenght_all_instances = len(all_lines)
    if lenght_all_instances > 0:
        all_lines.sort(key=lambda k: k[score_index], reverse=True)
    
    number_of_highest_ranked_instances_to_write = (lenght_all_instances*percentage_highest_scored_isntances)/100
    number_of_highest_ranked_instances_wrote = 0
    number_of_sig_instances_wrote = 0
    if lenght_all_instances>=limit_to_check:#only check when at least (limit_to_check) number of motif instances are found
        mean_scores = np.mean(all_scores)
        sd_scores = np.std(all_scores)#this helps including scores in case the scores follow a uniform distribution
        
        with open(bed_input_file, 'w') as bed_outfile:
            with open(sig_bed_output_file, 'w') as sig_bed_outfile:
                with open(sig_ranked_bed_output_file, 'w') as sig_ranked_bed_outfile:
                    for sl in all_lines:
                        #sl[1] = str(int(sl[1])+1) #in case of using --parse-genomic-coordinates and if the coordinates are generated from mergeBed one base might be needed to add to each motif 
                        #sl[2] = str(int(sl[2])+1)
                        line_to_write = '\t'.join(sl[2:5]) + '\t' + sl[0] + '\t' + str(sl[6]) + "\t"+ str(sl[7]) + '\t' + sl[5] + '\n' 
                        if sl[score_index] >= np.floor(mean_scores+(sd_scores*scores_sd_above_mean)):
                            sig_bed_outfile.write(line_to_write)
                            number_of_sig_instances_wrote+=1
                        if number_of_highest_ranked_instances_wrote<number_of_highest_ranked_instances_to_write:
                            sig_ranked_bed_outfile.write(line_to_write)
                            number_of_highest_ranked_instances_wrote+=1
                        bed_outfile.write(line_to_write)
                    
    print(sig_bed_output_file + "(len_all_instances,mean_scores,sd_scores,sig_instances,highest_ranked_instances):\t" + str(lenght_all_instances) + '\t' + str(mean_scores) + '\t' + str(sd_scores) + '\t' + str(number_of_sig_instances_wrote) + '\t' + str(number_of_highest_ranked_instances_wrote))
    return bed_input_file, sig_bed_output_file
    
        
def distribute_meme_pwms_to_single_pwm(jaspar_meme_pwms_input_file, ):
    '''This is a helper function that processes the given PWM file to create a PWM file for each TF model listed in the given combined-PWM file
    It retruns a list of the created PWM files'''
    print('generate single pwm')
    list_pwm_files = []
    pwms_lines = []
    pwm_out_dir = jaspar_meme_pwms_input_file+"_single"
    if not os.path.exists(pwm_out_dir):
        os.makedirs(pwm_out_dir)
    
    with open(jaspar_meme_pwms_input_file, 'r') as pwms_infile:
        pwms_lines =  pwms_infile.readlines()
    header_lines = ''.join(pwms_lines[0:9])
    lines_to_write = ""
    pwm_out_file_name = ""
    for l in pwms_lines[9::]:
        if l.startswith("MOTIF"):
            split_l = l.strip().split(' ')
            if lines_to_write!="" and pwm_out_file_name!="":
                #write_the_previous_pwm
                with open(pwm_out_file_name, 'w') as pwm_outfile:
                    pwm_outfile.write(lines_to_write)
                list_pwm_files.append(pwm_out_file_name)
            if len(split_l)==3:#From JASPAR2016 DB
                pwm_out_file_name = pwm_out_dir + '/' + split_l[2] +  "_" + split_l[1] + ".pwm"
                lines_to_write = header_lines + "\nMOTIF " +  split_l[2] +  "_" + split_l[1] + " " + split_l[2]+'\n'
            if len(split_l)==2: #From HOMERv10 DB
                pwm_out_file_name = pwm_out_dir + '/' + split_l[1].split('_')[0] +  "_" + split_l[1].split('_')[1].split('.')[2] + ".pwm"
                lines_to_write = header_lines + "\nMOTIF " + split_l[1].split('_')[0] +  "_" + split_l[1].split('_')[1].split('.')[2] + ' ' + split_l[1] + '\n'
        else:
            lines_to_write += l
    if lines_to_write!="" and pwm_out_file_name!="":#to write the last motif occurence
        #write_the_previous_pwm
        with open(pwm_out_file_name, 'w') as pwm_outfile:
            pwm_outfile.write(lines_to_write)
        list_pwm_files.append(pwm_out_file_name)
    #fimo --oc test_fimo --max-strand --verbosity 5 --parse-genomic-coord --text --skip-matched-sequence head_pwm.meme2 ../HOCOMOCOv10/CellInfo21Sep16AllENCODE_allChIPSeqDNase.fa > test_fimo.txt
    return list_pwm_files

def get_motif_instances_FIMO_singlepwm(pwm_file, 
                                       fasta_file, 
                                       instances_file_fimo_format, 
                                       instances_file_bed_format, 
                                       sig_instances_file_bed_format, 
                                       sig_ranked_bed_output_file, 
                                       threshold, 
                                       limit_to_check, 
                                       scores_sd_above_mean, 
                                       percentage_highest_scored_isntances):
    '''Runs fimo for a given PWM file and generate its motif instances
    Calls get_sig_motif_instances to generate the significant motif instances
    Returns a bed file containing significant motif instances'''
    if not os.path.exists(instances_file_fimo_format):
        fimo_stm = """fimo --max-strand --parse-genomic-coord --skip-matched-sequence --thresh "%s" "%s" "%s" > "%s" """%(str(threshold),pwm_file, fasta_file, instances_file_fimo_format)
        #awk 'BEGIN{FS=OFS="\t"}{if($6=="%s") print $0 >> "%s"}' %s""" %(cohort_value, cohort_file, mutations_input_file)
        os.system(fimo_stm)
        #fimo_out  = open(instances_file_fimo_format, 'w')
        #fimo_cmd = "fimo --max-strand --parse-genomic-coord --skip-matched-sequence --thresh" + " " + str(threshold) + " " + pwm_file + " " + fasta_file 
        
        #proc = subprocess.Popen(shlex.split(fimo_cmd), stdout=PIPE, stderr=STDOUT)
        #stdout, stderr = proc.communicate()
        #fimo_out.write(stdout)
        #fimo_out.close()
    if not os.path.exists(sig_instances_file_bed_format) or not os.path.exists(sig_ranked_bed_output_file):
        get_sig_motif_instances(instances_file_fimo_format, instances_file_bed_format, sig_instances_file_bed_format, sig_ranked_bed_output_file, limit_to_check, scores_sd_above_mean, percentage_highest_scored_isntances)
    
    return sig_instances_file_bed_format
    
def run_motifs_FIMO(pwm_files, 
                    fasta_file, 
                    output_dir, 
                    threshold, 
                    limit_to_check, 
                    scores_sd_above_mean, 
                    percentage_highest_scored_isntances, 
                    number_processes_to_run_in_parallel=2):
    '''Reads a given list of PWM files
    Calls get_motif_instances_FIMO_singlepwm to generate significant instances for each PWM file in parallel
    Returns a list of generated significant bed files'''
    generated_FIMO_files = []
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    pool = Pool(processes=number_processes_to_run_in_parallel)
    for pwm_file in pwm_files:
        motif_name = '.'.join(pwm_file.split('/')[-1].split('.')[0:-1])
        instances_file_fimo_format = output_dir+'/'+motif_name+".fimo"
        instances_file_bed_format = output_dir+'/'+motif_name+".bed"
        sig_instances_file_bed_format = output_dir+'/'+motif_name+"_sig.bed"
        sig_ranked_bed_output_file = output_dir+'/'+motif_name+"_highestranked" + str(percentage_highest_scored_isntances) + ".bed"
        if not os.path.exists(sig_instances_file_bed_format) or not os.path.exists(sig_ranked_bed_output_file):
            pool.apply_async(get_motif_instances_FIMO_singlepwm, 
                             args=(pwm_file, fasta_file, instances_file_fimo_format, instances_file_bed_format, sig_instances_file_bed_format, sig_ranked_bed_output_file, threshold, limit_to_check, scores_sd_above_mean, percentage_highest_scored_isntances), callback=generated_FIMO_files.append) 
    pool.close()
    pool.join()
    return generated_FIMO_files


def parse_args():
    '''Parse command line arguments'''
    print('parse')
    parser = argparse.ArgumentParser(description='Generate Motifs Fimo')
    parser.add_argument('--jaspar_meme_pwms_input_file', default='', help='')
    parser.add_argument('--genome_fasta_file', default='', help='')    
    parser.add_argument('--output_dir', default='', help='') 
    parser.add_argument('--pval_threshold', default=0.0001, help='',type=float)
    parser.add_argument('--limit_to_check', default=1, help='', type=int)   
    parser.add_argument('--scores_sd_above_mean', default=1.0, help='', type=float)   
    parser.add_argument('--percentage_highest_scored_isntances', default=10, help='',type=int) 
    parser.add_argument('--number_processes_to_run', default=16, help='', type=int) 
    
  

    

    
    return parser.parse_args(sys.argv[1:])
 

if __name__ == '__main__':
    
    args = parse_args()
       

    
    #print("Usage: python GenerateMotifsFIMO.py jaspar_meme_pwms_input_file.txt hg19.fa output_dir pval_threshold<float> limit_to_check<int> scores_sd_above_mean<float> percentage_highest_scored_isntances<int> #processes<int>")
  
    pwms_files = distribute_meme_pwms_to_single_pwm(jaspar_meme_pwms_input_file=args.jaspar_meme_pwms_input_file)
    run_motifs_FIMO(pwm_files=pwms_files, fasta_file=args.genome_fasta_file, output_dir=args.output_dir, threshold = float(args.pval_threshold), limit_to_check = int(args.limit_to_check), scores_sd_above_mean=float(args.scores_sd_above_mean), percentage_highest_scored_isntances = int(args.percentage_highest_scored_isntances), number_processes_to_run_in_parallel=int(args.number_processes_to_run))
    
    