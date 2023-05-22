# SnapHiC-D

## Identifying differential chromatin interactions from single cell Hi-C data

### Last edited: 05/2023. Version: 0.1.0

Find the preprint [here](https://www.biorxiv.org/content/10.1101/2022.08.05.502991v1).

SnapHiC-D is an extension of [SnapHiC](https://github.com/HuMingLab/SnapHiC) and requires SnapHiC's RWR step output for its input. For a faster version, use SnapHiC2, which is enabled by selecting "method="sliding_window".

### Requirements
SnapHiC-D was built using following Python packages.

1. Python 3.6.8
2. numpy 1.19.5
3. pandas 1.1.5
4. qnorm 0.8.1 (https://github.com/Maarten-vd-Sande/qnorm)
5. scipy 1.5.4
6. statsmodels 0.12.2
7. futures 3.0.5

### Running SnapHiC-D

The Shellscript is included in this github page. The required inputs variables are:

1. SnapHiC_D_dir : The directory of SnapHiC-D
2. group_A_dir : The directory of files for group A 
3. group_B_dir : The directory of files for group B
4. file_list_dir : The directory of file containing file name and its group information (i.e. A or B) with the 2 column header of "name group". One may use print_file_group_name.R code to print this file. The example of a file list is "OM_chr3_final.txt".
5. out_dir : The output directory
6. chr : chromosome number (i.e. chr3)
7. genome : genome (i.e. mm10 or hg19)
8. num_CPUs : The number of CPUs one would like to use. One can check how many CPUs are available by "lscpu". If num_CPUs = 1, the program will run as a  single processor. When using a HPC with job scheduler, make sure to ask for 1 node. 
9. bin_size : The resolution of bin size
10. fdr_threshold : FDR threshhold; the default value is 0.1
11. max_gap : The maximum distance gap; the default value is 101 (1MB)
12. min_gap : The minimum distance gap; the default value is 2 (2kb)

To run SnapHiC-D,
```
./run_SnapHiC_D.sh
```

### Contact Us
For any questions regarding this software, contact Ming Hu (hum@ccf.org) or Lindsay Lee (leeh7@ccf.org).
