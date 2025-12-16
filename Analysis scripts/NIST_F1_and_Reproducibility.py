#! /home/dsm_23110700129/miniconda3/envs/strelka_env/bin/python2.7
# -*- coding: UTF-8 -*-
import pandas as pd
import os
import itertools
import sys
import argparse

parser = argparse.ArgumentParser(description='Compare with NIST reference datasets and calculate F1 score')
parser.add_argument('-i', '--input_dir', help='Path to the input files, includes *.vcf.gz.')
parser.add_argument('-o', '--out_dir', help='Path to the output dir')
parser.add_argument('-s', '--ref_set',default='no', help='Set to yes if comparing with reference dataset, otherwise set to no')
parser.add_argument('-R','--region',required=False,help='Specify the region bed file, comparison will only be done within this region')
parser.add_argument('-t','--threds',required=False,default=4,help='Set the number of threads')
args = parser.parse_args()

input_dir = args.input_dir
ref_set = args.ref_set
out_dir = args.out_dir
region = args.region
thred = args.threds


reference = '/cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/Reference/GRCh38.d1.vd1.fa'

def unique_combinations(arr):
    combos = itertools.combinations(arr, 2)
    unique_combos = set(combos)
    return unique_combos

aim_arr = []

for i in os.listdir(input_dir):
    if i.endswith('.gz') or i.endswith('.vcf'):
        aim_arr.append(i)

print(aim_arr)
#out_dir = caller+'_'+year
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

HC_SNV = '/cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/HCC1395_DNA/NIST_high_confidence/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz'

if ref_set == 'yes':
    for a in aim_arr:
        #In high-confidence region
        name_SNV = a.replace('.vcf','').replace('.gz','')+'.SNV_HC'
        name_Indel = a.replace('.vcf.gz','')+'.Indel_HC'
        #SNV
        os.system('hap.py {0} {1} -r {2} -T {3} -o {4} --threads {5}'.
                  format(HC_SNV,
                         input_dir+'/'+a,
                         reference,
                         region,
                         out_dir+'/'+name_SNV,thred))
        

else:
    for out in unique_combinations(aim_arr):
        a,b = out[0],out[1]
        name = a.replace('.vcf.gz','')+'VS'+b.replace('.vcf.gz','')
        print('-------Processing '+name+'------')
        #In the whole genome
        if not region:
            print('Whole genome........')
            os.system('hap.py {0} {1} -r {2} -o {3} --threads {4}'.format(
                input_dir+'/'+a,
                input_dir+'/'+b,
                reference,
                out_dir+'/'+name,thred))
        else:
            print('Target region..........')
        #Only in high-confidence region
            os.system('hap.py {0} {1} -r {2} -R {3} -o {4} --threads {5}'.format(
                input_dir+'/'+a,
                input_dir+'/'+b,
                reference,
                region,
                out_dir+'/'+name+'_HCR',thred))


file_list=[]
for i in os.listdir(out_dir):
    if i.endswith('.summary.csv'):
        stats = pd.read_csv(out_dir+'/'+i,sep=',')
        name = i.split('.summary')[0]
        stats_sub = stats[stats['Filter']=='ALL']
        stats_sub['source']=name
        print(stats_sub)

        file_list.append(stats_sub)


stats_df = pd.concat(file_list)
stats_df.to_csv(out_dir+'/'+'F1score_stats_total_NIST.csv',index=None)
os.system('rm {}/*.vcf*'.format(out_dir))
os.system('rm {}/*.csv.gz'.format(out_dir))
