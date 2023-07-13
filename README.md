# SnapHiC-D

## Identifying differential chromatin interactions from single cell Hi-C data

Find the preprint [here](https://www.biorxiv.org/content/10.1101/2022.08.05.502991v1).

SnapHiC-D is an extension of [SnapHiC](https://github.com/HuMingLab/SnapHiC) and requires SnapHiC's RWR step output for its input. For a faster version, use SnapHiC2, which is enabled by selecting "method="sliding_window".

### Install SnapHiC-D

Install SnapHiC-D through `pip`:

```
conda create --name SnapHiC_D_env python==3.6.8
conda activate SnapHiC_D_env
pip install SnapHiC-D
```

### Requirements
SnapHiC-D was built using following Python packages.

1. Python 3.6.8
2. numpy 1.19.0
3. pandas 1.1.5
4. qnorm 0.8.1 (https://github.com/Maarten-vd-Sande/qnorm)
5. scipy 1.5.4
6. statsmodels 0.12.2
7. futures 3.0.5
8. click 7.1.2

### Running SnapHiC-D

Activate the python environment with SnapHiC-D installed and enter the following in the terminal:

```
SnapHiC-D diff-loops -i group_A_dir -j group_B_dir -o out_dir -c chr -n num_CPUs\
                     -b genome_region_path -g genome_transcript_path\
                     --binsize bin_size --fdr_threshold fdr_threshold\
                     --mini_gap min_gap --maxi_gap max_gap
```

The required inputs variables are:

1. group_A_dir : The directory of files for group A 
2. group_B_dir : The directory of files for group B
3. out_dir : The output directory
4. chr : chromosome number (i.e. chr3)
5. num_CPUs : The number of CPUs one would like to use. One can check how many CPUs are available by "lscpu". If num_CPUs = 1, the program will run as a  single processor. When using a HPC with job scheduler, make sure to ask for 1 node. 
6. genome_region_path: the path of *mm10_filter_regions.txt* or *hg19_filter_regions.txt*, depending on the reference genome. These files are provided in the *ext* folder.
7. genome_transcript_path: the path of *mm10.refGene.transcript.TSS.061421.txt* or *hg19.refGene.transcript.TSS.061421.txt*, depending on the reference genome. These files are provided in the *ext* folder.
8. bin_size : The resolution of bin size
9. fdr_threshold : FDR threshhold; the default value is 0.1
10. min_gap : The minimum distance gap; the default value is 2 (2kb)
11. max_gap : The maximum distance gap; the default value is 101 (1MB)

We have provided input example data of 94 mouse embryonic stem cells (mESC) and 188 mouse neuron progenitor cells (NPCs) in zipped folders to test SnapHiC-D. These are the trimmed RWR results from SnapHiC around the 200Kb region of *Sox2* locus - chr3:34,601,000–34,806,000 (ref: mm10). To run SnapHiC-D, type

```
SnapHiC-D diff-loops -i group_A_dir -j group_B_dir -o output -c chr3 -n 2\
                     -b "ext/mm10_filter_regions.txt"\
                     -g "ext/mm10.refGene.transcript.TSS.061421.txt"
```

A directory named *output* will be created with the following files inside:

1. *output/combined_results_chr3.txt*: *T*-test results of bin pairs.
2. *output/DI_FDR0.1_T2_Test_chr3.txt*: filtered results based on FDR and the *T* statistic.

### Contact Us
For any questions regarding this software, contact Ming Hu (hum@ccf.org), Lindsay Lee (leeh7@ccf.org), or Hongyu Yu (hongyuyu@unc.edu).
