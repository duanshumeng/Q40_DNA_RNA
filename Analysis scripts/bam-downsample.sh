#1. Count the total number of bases in the exon sequencing interval
zcat hg38_exome_comp_spikein_v2.0.2_targets_sorted.re_annotated_noy_cut.bed.gz | awk 'BEGIN { totalLength = 0 } { totalLength += $3 - $2 } END { print "Total Length: " totalLength }'
Total Length: 37352797

#2. Count the number of bam files
total_num=`samtools view -@ 10 -c HCC1395BL_Fresh_1.sorted.deduped.recaled.bam`

#3. Calculate the percentage of reads at the specified depth in the total
freq30X=`awk -v vv=$total_num 'BEGIN {print 7470559/(vv/2)}'`
freq60X=`awk -v vv=$total_num 'BEGIN {print 14941119/(vv/2)}'`
freq90X=`awk -v vv=$total_num 'BEGIN {print 22411678/(vv/2)}'`
freq120X=`awk -v vv=$total_num 'BEGIN {print 29882237/(vv/2)}'`

#4. Only Germline samples are processed
freq10X=`awk -v vv=$total_num 'BEGIN {print 4980373/(vv/2)}'`
freq15X=`awk -v vv=$total_num 'BEGIN {print 7470560/(vv/2)}'`
freq20X=`awk -v vv=$total_num 'BEGIN {print 9960746/(vv/2)}'`
freq25X=`awk -v vv=$total_num 'BEGIN {print 12450932/(vv/2)}'`

export PATH="/cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/miniconda3/bin:$PATH"
source activate mapping_V2
cd /cpfs01/projects-HDD/cfff-e44ef5cf7aa5_HDD/dsm_23110700129/HCC1395_DNA/Elemenet.vs.Illumina/MGI_Q40/MGI_WES/Bam_downsample
input_file="wes_downsample_90X.txt"
mapfile -t samples < <(awk -F'\t' -v OFS=',' 'NR > 1 {print $1 "," $2 "," $3 "," $4 "," $5}' "$input_file")  
for sample_info in "${samples[@]}"; do                                                       
        IFS=',' read -ra sample_fields <<< "$sample_info"                               
        sample_id="${sample_fields[0]}"    
        bam="${sample_fields[1]}"                                                           
        depth="${sample_fields[2]}"  
        n="${sample_fields[3]}"  
        out_bam="${sample_fields[4]}" 
        echo $sample_id
		set -o pipefail
        set -e
        nt=8
        (total_num=`samtools view -@ $nt -c ${bam}`
        echo $total_num
        num=${n}
        #freq=$(awk -v vv="$total_num" -v aa="$num" 'BEGIN {print aa / vv}')
        seed_freq=$(awk -v vv="$total_num" -v aa="$num" 'BEGIN {print 100+(aa / vv)}')
        samtools view -@ $nt -s $seed_freq $bam -b -o $out_bam
        samtools index -b $out_bam) &
 done  
 wait 
 
 echo "Downsampling FInished...."