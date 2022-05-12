# SnapHiC-D

## Identifying differential chromatin interactions from single cell Hi-C data

### Last edited: 05/2022. Version: 0.1.0


SnapHiC-D is an extension of SnapHiC (https://github.com/HuMingLab/SnapHiC) and requires SnapHiC's RWR step output for its input. For a faster version, use SnapHiC2, which is enabled by selecting "method="sliding_window".

### Requirements
SnapHiC-D was built using following Python packages.

1. Python 3.6.8
2. numpy 1.19.5
3. pandas 1.1.5
4. qnorm 0.8.1
5. scipy 1.5.4
6. statsmodels 0.12.2
7. futures 3.0.5

### Running SnapHiC-D

The shellscript is included in the github page. The required inputs variables are:

1. group_A_dir : The directory of files for group A 
2. group_B_dir : The directory of files for group B
3. file_list_dir : The directory of file containing file name and its group information (i.e. A or B). One may use print_file_group_name.R code to print this file 
4. out_dir : The output directory
5. chr : chromosome number (i.e. chr3)
6. num_CPUs : The number of CPUs one would like to use. One can check how many CPUs are available by "lscpu". If num_CPUs = 1, the program will run as single processor. 
7. bin_size : The resolution of bin size
8. fdr_threshold : FDR threshhold; the default value is 0.1
9. max_gap : The maximum distance gap; the default value is 101 (1MB)
10. min_gap : The minimum distance gap; the default value is 2 (2kb)
