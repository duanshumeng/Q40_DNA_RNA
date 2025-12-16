#Extract VAF information
python /cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/HCC1395_DNA/python_scripts/calculate_VAF.py -i WES_HCC1395_1_ELE_10G.90X.TNseq.vcf -c TNseq -R /cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/HCC1395_DNA/Validation_by_WES_element/hg38_exome_comp_spikein_v2.0.2_targets_sorted.re_annotated_noy_cut.HCR.bed

#Combine VAF information,get mena and sd
def combine_VAF(input_dir,data_type):
    vafs = subprocess.check_output("ls {0}/*{1}*.VAF.csv".format(input_dir,data_type),shell=True)
    vaf_files = vafs.decode('utf-8').split('\n')[:-1]
    vaf_df = pd.read_csv(vaf_files[0])
    vaf_df = vaf_df[['VAF', 'tag', 'Type']]
    smp = os.path.basename(vaf_files[0]).split('.')[0]
    vaf_df.columns = [smp+' VAF', 'tag', smp+' Type']
    for i in vaf_files[1:]:
        print(i)
        df = pd.read_csv(i)
        df = df[['VAF', 'tag', 'Type']]
        smp = os.path.basename(i).split('.')[0]
        df.columns = [smp+' VAF', 'tag', smp+' Type']
        vaf_df=pd.merge(vaf_df,df,on='tag', how='inner')
    for i in vaf_df.columns:
        if 'VAF' in i:
            vaf_df[i] = vaf_df[i].apply(lambda x:float(x.split(',')[0]))

    vaf_df['VAF Mean'] = vaf_df[[i for i in vaf_df.columns if 'VAF' in i]].mean(axis=1)
    vaf_df['VAF STD'] = vaf_df[[i for i in vaf_df.columns if 'VAF' in i]].std(axis=1)
    vaf_df['CV'] = vaf_df['VAF STD']/vaf_df['VAF Mean']
    vaf_df.to_csv('VAF.{}.csv'.format(data_type),index=None)

input_dir='/cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/HCC1395_DNA/Elemenet.vs.Illumina/Data_Performance/HCC1395_WES/TNseq_120X'
combine_VAF(input_dir,'ELE')
combine_VAF(input_dir,'ILM')
