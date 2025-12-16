from __future__ import division
import pandas as pd
import sys, argparse, os
import fileinput
import re

# input arguments
parser = argparse.ArgumentParser(description="this script is to extract mendelian concordance information")

parser.add_argument('-LCL5', '--LCL5', type=str, help='LCL5 family info',  required=True)
parser.add_argument('-LCL6', '--LCL6', type=str, help='LCL6 family info',  required=True)
parser.add_argument('-genotype', '--genotype', type=str, help='Genotype information of a set of four family members',  required=True)
parser.add_argument('-family', '--family', type=str, help='family name',  required=True)


args = parser.parse_args()
lcl5 = args.LCL5
lcl6 = args.LCL6
genotype = args.genotype
family = args.family


# output file
family_name = family + '.txt'

family_file = open(family_name,'w')

summary_name = family + '.summary.txt'

summary_file = open(summary_name,'w')

# input files
lcl5_dat = pd.read_table(lcl5)
lcl6_dat = pd.read_table(lcl6)
genotype_dat = pd.read_table(genotype)
merged_df = pd.merge(lcl5_dat, lcl6_dat,  how='outer', left_on=['#CHROM','POS'], right_on = ['#CHROM','POS'])
merged_genotype_df = pd.merge(merged_df, genotype_dat,  how='outer', left_on=['#CHROM','POS'], right_on = ['#CHROM','POS'])
merged_genotype_df = merged_genotype_df.dropna()

merged_genotype_df_sub = merged_genotype_df.iloc[:,[0,1,23,24,29,30,31,32,7,17]]
merged_genotype_df_sub.columns = ['CHROM', 'POS', 'REF', 'ALT','LCL5','LCL6','LCL7','LCL8', 'TRIO5', 'TRIO6']

indel_sister_same = 0
indel_sister_diff = 0
indel_family_all = 0
indel_family_mendelian = 0

snv_sister_same = 0
snv_sister_diff = 0
snv_family_all = 0
snv_family_mendelian = 0


for row in merged_genotype_df_sub.itertuples():
	# indel or snv
	if ',' in row.ALT:
		alt = row.ALT.split(',')
		alt_len = [len(i) for i in alt]
		alt_max = max(alt_len)
	else:
		alt_max = len(row.ALT)
	alt = alt_max
	ref = row.REF
	if len(ref) == 1 and alt == 1:
		cate = 'SNV'
	elif len(ref) > alt:
		cate = 'INDEL'
	elif alt > len(ref):
		cate = 'INDEL'
	elif len(ref) == alt:
		cate = 'INDEL'
	# sister
	if row.LCL5 == row.LCL6:
		if row.LCL5 == './.':
			mendelian = 'noInfo'
			sister_count = "no"
		elif row.LCL5 == '0/0':
			mendelian = 'Ref'
			sister_count = "no"
		else:
			mendelian = '1'
			sister_count = "yes_same"
	else:
		mendelian = '0'
		if (row.LCL5 == './.' or row.LCL5 == '0/0') and (row.LCL6 == './.' or row.LCL6 == '0/0'):
			sister_count = "no"
		else:
			sister_count = "yes_diff"	
	# family trio5
	if row.LCL5 == row.LCL7 == row.LCL8 == './.':
		mendelian = mendelian + ':noInfo'
	elif row.LCL5 == row.LCL7 == row.LCL8 == '0/0':
		mendelian = mendelian + ':Ref'
	elif pd.isnull(row.TRIO5) == True:
		mendelian = mendelian + ':unVBT'
	else:
		mendelian = mendelian + ':' + row.TRIO5.split('=')[1]
	# family trio6
	if row.LCL6 == row.LCL7 == row.LCL8 == './.':
		mendelian = mendelian + ':noInfo'
	elif row.LCL6 == row.LCL7 == row.LCL8 == '0/0':
		mendelian = mendelian + ':Ref'
	elif pd.isnull(row.TRIO6) == True:
		mendelian = mendelian + ':unVBT'
	else:
		mendelian =  mendelian + ':' + row.TRIO6.split('=')[1]
	# not count into family
	if (row.LCL5 == './.' or row.LCL5 == '0/0') and (row.LCL6 == './.' or row.LCL6 == '0/0') and (row.LCL7 == './.' or row.LCL7 == '0/0') and (row.LCL8 == './.' or row.LCL8 == '0/0'):
		mendelian_count = "no"
	else:
		mendelian_count = "yes"
	outline = str(row.CHROM) + '\t' + str(row.POS) + '\t' + str(row.REF) + '\t' + str(row.ALT) + '\t' + str(cate) + '\t' + str(row.LCL5) + '\t' + str(row.LCL6) + '\t' + str(row.LCL7) + '\t' + str(row.LCL8) + '\t' + str(row.TRIO5) + '\t' + str(row.TRIO6) + '\t' + str(mendelian) + '\t' + str(mendelian_count) + '\t' + str(sister_count) + '\n'
	family_file.write(outline)
	if cate == 'SNV':
		if sister_count == 'yes_same':
			snv_sister_same += 1
		elif sister_count == 'yes_diff':
			snv_sister_diff += 1
		else:
			pass
		if mendelian_count == 'yes':
			snv_family_all += 1
		else:
			pass
		if mendelian == '1:1:1':
			snv_family_mendelian += 1
		elif mendelian == 'Ref:1:1':
			snv_family_mendelian += 1
		else:
			pass		
	elif cate == 'INDEL':
		if sister_count == 'yes_same':
			indel_sister_same += 1
		elif sister_count == 'yes_diff':
			indel_sister_diff += 1
		else:
			pass	
		if mendelian_count == 'yes':
			indel_family_all += 1
		else:
			pass
		if mendelian == '1:1:1':
			indel_family_mendelian += 1
		elif mendelian == 'Ref:1:1':
			indel_family_mendelian += 1
		else:
			pass

snv_sister = snv_sister_same/(snv_sister_same + snv_sister_diff)
indel_sister = indel_sister_same/(indel_sister_same + indel_sister_diff)
snv_quartet = snv_family_mendelian/snv_family_all
indel_quartet = indel_family_mendelian/indel_family_all
outcolumn =  'Family\tTotal_Variants\tMendelian_Concordant_Variants\tMendelian_Concordance_Rate\n'
indel_outResult = family + '.INDEL' + '\t' + str(indel_family_all) + '\t' + str(indel_family_mendelian) + '\t' + str(indel_quartet) + '\n'
snv_outResult = family + '.SNV' + '\t' + str(snv_family_all) + '\t' + str(snv_family_mendelian)  + '\t' + str(snv_quartet) + '\n'
summary_file.write(outcolumn)
summary_file.write(indel_outResult)
summary_file.write(snv_outResult)



