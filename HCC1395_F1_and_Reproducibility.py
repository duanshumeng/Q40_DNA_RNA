#! /home/dsm_23110700129/miniconda3/envs/strelka_env/bin/python2.7
# -*- coding: UTF-8 -*-
import pandas as pd
import os
import itertools
import sys
import argparse

parser = argparse.ArgumentParser(description='Compare with SEQC2 reference datasets and calculate F1 score')
parser.add_argument('-i', '--input_dir', help='Path to the input files, includes *.vcf.gz.')
parser.add_argument('-o', '--out_dir', help='Path to the output dir')
parser.add_argument('-s', '--ref_set',default='no', help='如果是与参考数据集比较，则设置为yes，否则为no')
parser.add_argument('-R','--region',required=False,help='指定区间的bed文件，仅在指定区间内进行比较')
args = parser.parse_args()

input_dir = args.input_dir
ref_set = args.ref_set
out_dir = args.out_dir
region = args.region
#如果是和seqc2的高置信数据集比较，则设置为yes，否则为no
#command demo

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

#SEQC2:
#HC_SNV ="/cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/HCC1395_DNA/SEQC2_high_confidence/high-confidence_in_HC_regions_v1.2.1_sort_SNV_Indels.vcf.gz"
HC_SNV = '/cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/HCC1395_DNA/SEQC2_high_confidence/high-confidence_in_HC_regions_v1.2.1_lowVAF.vcf.gz'
#region ="/cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/HCC1395_DNA/SEQC2_high_confidence/High-Confidence_Regions_v1.2.sorted.bed"
#SEQC2:
#HC_SNV = "/cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/HCC1395_DNA/High_confidence_datasets/high-confidence_dataset2_V2/high-confidence_sSNV_sIndel_v2.HCR.clean.sort.vcf.gz"
#HC_SNV = '/cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/HCC1395_DNA/High_confidence_datasets/high-confidence_dataset2_V2/high-confidence_sSNV_sIndel_v2.HCR.lowVAF.vcf.gz'
#region = "/cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/HCC1395_DNA/High_confidence_datasets/high-confidence_dataset2_V2/PGx_High-Confidence_Regions_v1.6.bed"

if ref_set == 'yes':
    for a in aim_arr:
        #仅在高置信区间内比较
        name_SNV = a.replace('.vcf','').replace('.gz','')+'.SNV_HC'
        name_Indel = a.replace('.vcf.gz','')+'.Indel_HC'
        #Indel
        os.system('som.py {0} {1} -r {2} -R {3} -o {4}'.
                  format(HC_SNV,
                         input_dir+'/'+a,
                         reference,
                         region,
                         out_dir+'/'+name_SNV))
        # # Indel
        # os.system('som.py {0} {1} -r {2} -R {3} -o {4}'.
        #           format(HC_Indel,
        #                  input_dir+'/'+a,
        #                  reference,
        #                  region,
        #                  out_dir+'/'+name_Indel))

        

else:
    for out in unique_combinations(aim_arr):
        a,b = out[0],out[1]
        name = a.replace('.vcf.gz','')+'VS'+b.replace('.vcf.gz','')
        print('-------Processing '+name+'------')
        #全基因组范围内比较
        if not region:
            print('Whole genome........')
            os.system('som.py {0} {1} -r {2} -o {3}'.format(
                input_dir+'/'+a,
                input_dir+'/'+b,
                reference,
                out_dir+'/'+name))
        else:
            print('Target region..........')
        #仅在高置信区间内比较
            os.system('som.py {0} {1} -r {2} -R {3} -o {4}'.format(
                input_dir+'/'+a,
                input_dir+'/'+b,
                reference,
                region,
                out_dir+'/'+name+'_HCR'))


file_list=[]
for i in os.listdir(out_dir):
    if i.endswith('.stats.csv'):
        stats = pd.read_csv(out_dir+'/'+i,sep=',')
        name = i.split('.stats')[0]
        stats_sub = stats[['type','total.truth','total.query','tp','fp','fn','recall','precision']]
        stats_sub['source']=name
        stats_sub['F1 score'] = stats_sub.apply(lambda df:2 * (df['precision'] * df['recall']) / (df['precision'] + df['recall']) if (df['precision'] + df['recall']) != 0 else 0,axis=1)
        stats_sub['Reproducibility'] = stats_sub.apply(lambda df:(float(df['tp'])/df['total.truth'] + float(df['tp'])/df['total.query'])/2 if (df['precision'] + df['recall']) != 0 else 0,axis=1)
        file_list.append(stats_sub)

stats_df = pd.concat(file_list)
stats_df.to_csv(out_dir+'/'+'F1score_stats_total.csv',index=None)
os.system('rm {}/*.txt'.format(out_dir))
os.system('rm {}/*.vcf'.format(out_dir))
os.system('rm {}/*.stats'.format(out_dir))
