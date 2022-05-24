#!/bin/bash

############################################################################
###                            User Variables                            ###
############################################################################

SnapHiC_D_dir="/home/leeh7/snapHiC-diff/python"
group_A_dir="/home/leeh7/RWR/ODC_RWR/ODC_rwr_chr3/"
group_B_dir="/home/leeh7/RWR/MG_final_RWR/MG_rwr_chr3/"
file_list_dir="/home/leeh7/snapHiC-diff/file_lists/ODC_MG/OM_chr3_final.txt"
out_dir="/home/leeh7/snapHiC-diff/python/output"
chr="chr3"
genome="hg19"
num_CPUs=5
bin_size=10000
fdr_threshold=0.1
max_gap=6
min_gap=2



############################################################################
if [ ! -d "$out_dir" ]; then
  mkdir $out_dir
fi

if [ ! -d "${out_dir}/tempfile" ]; then
  mkdir ${out_dir}/tempfile
fi

python ${SnapHiC_D_dir}/snapHiC_diff.py -s $SnapHiC_D_dir -i $group_A_dir -j $group_B_dir -a $file_list_dir -o $out_dir -c $chr -g $genome -n $num_CPUs --binsize $bin_size --fdr_threshold $fdr_threshold --maxi_gap $max_gap --mini_gap $min_gap
#python test_pars.py -i $group_A_dir -j $group_B_dir -a $file_list_dir -o $out_dir -c $chr -n $num_CPUs --binsize $bin_size --fdr_threshold $fdr_threshold --maxi_gap $max_gap --mini_gap $min_gap


if [ "$(ls -A ${out_dir}/tempfile)" ]; then
     cd ${out_dir}/tempfile; awk 'FNR==1 && NR!=1{next;}{print}' DI_FDR10_T2_Test_${chr}_GAP*.txt > ../DI_FDR10_T2_Test_all${chr}.txt; awk 'FNR==1 && NR!=1{next;}{print}' all_results_${chr}_GAP*.txt > ../combined_all${chr}_results.txt
 
else
    echo Program did not run successfully
fi