import numpy as np
import pandas as pd
import qnorm
from scipy import stats
import statsmodels.stats.multitest as multi
import concurrent.futures
from functools import partial
import argparse

import sys
import os



def main():
    parser = create_parser()
    args = parser.parse_args()
    
    CHR =  args.chrom
    genome = args.gen
    SnapHiC_dir = args.indirSHD
    dir_in_A = args.indirA
    dir_in_B = args.indirB
    ann_dir = args.anndir
    
    out_dir = args.outdir
    
    min_gap = args.mini_gap
    max_gap = args.maxi_gap
    
    BINSIZE = args.binsize
    
    fdr_t = args.fdr_threshold
    
    max_worker = args.num_cpus
    
    b_ID, gID = get_ID(CHR,BINSIZE,SnapHiC_dir,genome)
    GAP=list(range(min_gap, max_gap))
    out = []
    final_out = []
    with concurrent.futures.ProcessPoolExecutor(max_worker) as executor:
        for result in executor.map(partial(rest_of_steps, chr_list=CHR, dir_in_B=dir_in_B, dir_in_A=dir_in_A, \
                                            ann_dir=ann_dir, out_dir=out_dir, b_ID=b_ID, gID=gID, BINSIZE=BINSIZE, fdr_t=fdr_t) , GAP):
            
            print("..", flush=True) 
            
    print("Exiting program", flush=True) 

def g_filter(x):
    counts = np.count_nonzero(x > 1.96, axis=1)
    return counts


def get_usize(GAP, ann_dir, dir_in_A, dir_in_B,b_ID,gID,BINSIZE):
    ann = pd.read_csv(ann_dir, header='infer', sep=' ', lineterminator='\n')
    if(ann['group'][1] == "A"):
        x = pd.read_csv(dir_in_A + ann['name'][1], header=None, sep="\t")
    else:
        x = pd.read_csv(dir_in_B + ann['name'][1], header=None, sep="\t")

    x.columns = ["chrnum", "from", "from.1", "chrnum.1", "to", "to.1", "value"]
    x[["from", "from.1", "to", "to.1"]] = x[["from", "from.1", "to", "to.1"]]/BINSIZE

    y = x[x['to']-x['from'] == GAP]
    z = y[(~y['from'].isin(b_ID)) & (~y['to'].isin(b_ID))]
    u = z[(z['from'].isin(gID)) | (z['to'].isin(gID))]

    u_size = u.shape[0]
    return u_size


def get_ID(CHR,BINSIZE,SnapHiC_dir,genome):
    if(genome == "hg19"):
        b = pd.read_csv(os.path.join(SnapHiC_dir, "ext/hg19_filter_regions.txt"), sep="\t", header=None)
        g0 = pd.read_csv(os.path.join(SnapHiC_dir, "ext/hg19.refGene.transcript.TSS.061421.txt"), sep="\t", header='infer')
    if(genome == "mm10"):
        b = pd.read_csv(os.path.join(SnapHiC_dir, "ext/mm10_filter_regions.txt"), sep="\t", header=None)
        g0 = pd.read_csv(os.path.join(SnapHiC_dir, "ext/mm10.refGene.transcript.TSS.061421.txt"), sep="\t", header='infer')      
      
    b.columns = ['chr_name', 'x0', 'x1', 'x2', 'x3', 'x4', 'x5', 'x6']
    b1 = b.loc[b['chr_name'] == CHR]
    b_ID = b1.iloc[:, 2]/BINSIZE

    
    g = g0[g0['chr'] == CHR]
    g['TSS.Bin'] = np.ceil(g['TSS']/BINSIZE)
    gID = np.unique(g['TSS.Bin'])
    return b_ID, gID


def rest_of_steps(GAP, chr_list, dir_in_B, dir_in_A, ann_dir, out_dir, b_ID, gID, BINSIZE,fdr_t):
    print("Starting with GAP=", GAP)     
    CHR=chr_list
    ann = pd.read_csv(ann_dir, header='infer', sep=' ', lineterminator='\n')
    num_A = sum(ann['group'] == "A")
    num_B = sum(ann['group'] == "B")

    u_size = get_usize(GAP, ann_dir, dir_in_A, dir_in_B,b_ID,gID,BINSIZE)
    rec = np.zeros((u_size, len(ann)))
    
    
    for ID in range(0, len(ann)):
        if(ann['group'][ID] == "A"):
            x = pd.read_csv(dir_in_A + ann['name'][ID], header=None, sep="\t")
        else:
            x = pd.read_csv(dir_in_B + ann['name'][ID], header=None, sep="\t")

        x.columns = ["chrnum", "from", "from.1",
            "chrnum.1", "to", "to.1", "value"]
        x[["from", "from.1", "to", "to.1"]] = x[[
            "from", "from.1", "to", "to.1"]]/BINSIZE

        y = x[x['to']-x['from'] == GAP]
        z = y[(~y['from'].isin(b_ID)) & (~y['to'].isin(b_ID))]
        u = z[(z['from'].isin(gID)) | (z['to'].isin(gID))]

        rec[:, ID] = u['value']
        

    recA = rec[:, np.where(ann['group'] == "A")[0]]
    recB = rec[:, np.where(ann['group'] == "B")[0]]

    countA = g_filter(recA) > num_A*.1
    countB = g_filter(recB) > num_B*.1

    index = np.array([a or b for a, b in zip(countA, countB)]).astype(int)
    u_combine = np.c_[u[['from','from.1','to','to.1']], rec]
    final = u_combine[index == 1, ]
    
    if final.size == 0:
        return
    
    rec_trim = np.delete(final, np.s_[0:4], axis=1)

    rec_df = pd.DataFrame(rec_trim)
    rec_qq = qnorm.quantile_normalize(rec_df, axis=1)
    rec_qq_np = rec_qq.to_numpy()

    rec_qq_A = rec_qq_np[:, np.where(ann['group'] == "A")[0]]
    rec_qq_B = rec_qq_np[:, np.where(ann['group'] == "B")[0]]

    out = pd.DataFrame(columns=['x1','x2', 'y1','y2', 'mean.A',
                       'mean.B', 'Tstat', 'Ttest.Pvalue'])

    out['x1'], out['x2'] = final[:, 0], final[:, 1]
    out['y1'],out['y2'] = final[:,2], final[:,3]
    if all(np.var(rec_qq_A, ddof=1, axis=0) > 0) is True and all(np.var(rec_qq_B, ddof=1, axis=0) > 0) is True:
        t_stat, p_val = stats.ttest_ind(
            rec_qq_A, rec_qq_B, axis=1, equal_var=False)
        out['mean.A'] = np.mean(rec_qq_A, axis=1)
        out['mean.B'] = np.mean(rec_qq_B, axis=1)
        out['Tstat'] = t_stat
        out['Ttest.Pvalue'] = p_val
        rej, pval_corr = multi.multipletests(
            out['Ttest.Pvalue'], method='fdr_bh')[:2]
        out['fdr'] = pval_corr
#        out['x1'],out['x2'],out['y1'],out['y2'] = out['x1']-1,out['x2']-1,out['y1']-1,out['y2']-1
        out['x1'],out['x2'],out['y1'],out['y2'] = out['x1']*BINSIZE,out['x2']*BINSIZE,out['y1']*BINSIZE,out['y2']*BINSIZE
        out['d'] =out['y1']-out['x1']
        out.insert(0, "chr", CHR)
        final_out = out[(out['fdr'] < fdr_t) & (abs(out['Tstat']) > 2)]
        
      
        newfilename1 = os.path.join(out_dir, "tempfile/all_results_%s_GAP%s.txt" % (CHR,GAP)) 
        newfilename2 = os.path.join(out_dir, "tempfile/DI_FDR10_T2_Test_%s_GAP%s.txt" % (CHR,GAP)) 
        out.to_csv(newfilename1,
                        header=True, index=None, sep='\t')
        final_out.to_csv(newfilename2, header=True, index=None, sep='\t')
        
    else:
        print("Variance is not >0")

   # print("Done with GAP=", GAP, flush=True)        
    return out, final_out

def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--indirSHD', action = 'store', required = True, \
                        help = 'SnapHiC-D directory')
    parser.add_argument('-i', '--indirA', action = 'store', required = True, \
                        help = 'input group A directory')
    parser.add_argument('-j', '--indirB', action = 'store', required = True, \
                        help = 'input group B directory')                    
    parser.add_argument('-o', '--outdir', action = 'store', \
                        required = True, help = 'output directory')
    parser.add_argument('-a', '--anndir', action = 'store', \
                        required = True, help = 'file list directory')                    
    parser.add_argument('-c', '--chrom', action = 'store', help = 'chromosome to process', \
                        required = True)
    parser.add_argument('-g', '--gen', action = 'store', help = 'genome - mm10 or hg19', \
                        required = True)                        
    parser.add_argument('--binsize', type = int, help = 'bin size used for binning the reads', \
                        required = False, default = 1e4)
    parser.add_argument('-n', '--num_cpus', type = int , action = 'store', required = True, \
                        help = 'number of CPUS used for parallel computing')
    parser.add_argument('--fdr_threshold', default = 0.1, type = float, required = False, \
                        help = 'FDR threshold used for candidate peak detection')
    parser.add_argument('--maxi_gap', type = int, help = 'maximum gap between the read pairs', \
                        default = 100, required = False)
    parser.add_argument('--mini_gap', type = int, help = 'minimum gap between the read pairs', \
                        default = 2, required = False)
              
    return parser
    


    
if __name__ == "__main__":
    main()
