import os
import pandas as pd
import numpy as np
from scipy import stats
from qnorm import quantile_normalize
from statsmodels.stats import multitest as multi
from concurrent.futures import ProcessPoolExecutor
from functools import partial


def get_genome_ids(bid_path, gid_path, BINSIZE, chr_id):
    b = pd.read_csv(bid_path, header=None, sep="\t")
    bids = pd.unique((b[2]/BINSIZE)[b[0]==chr_id].astype("int"))
    g = pd.read_csv(gid_path, sep="\t")
    gids = np.ceil(g['TSS']/BINSIZE)[g["chr"]==chr_id]
    gids = pd.unique(gids).astype("int")
    return bids, gids


def read_one_rwr_file(rwr_path, BINSIZE, bids, gids, GAP, full=False):
    rwr_cols = ["chr1", "start1", "start2", "chr2", "end1", "end2", "val"]
    loc_cols = ["start1", "start2", "end1", "end2"]
    rwr_result = pd.read_csv(rwr_path, names=rwr_cols, sep="\t")
    rwr_result[loc_cols] = rwr_result[loc_cols]/BINSIZE
    rwr_subset = rwr_result[
        (rwr_result["end1"] - rwr_result["start1"] == GAP) & 
        ((~rwr_result["start1"].isin(bids))&(~rwr_result["end1"].isin(bids))) &
        (rwr_result["start1"].isin(gids) | rwr_result["end1"].isin(gids))
    ]
    return rwr_subset[loc_cols] if full else rwr_subset["val"].values


def single_gap_size(gap, chr_id, dire1, dire2, bids, gids, out_dire, BINSIZE):
    rwr_paths1 = [os.path.join(dire1, f) for f in os.listdir(dire1)]
    rec1 = np.stack([read_one_rwr_file(f, BINSIZE, bids, gids, gap) for f in rwr_paths1])
    rwr_paths2 = [os.path.join(dire2, f) for f in os.listdir(dire2)]
    rec2 = np.stack([read_one_rwr_file(f, BINSIZE, bids, gids, gap) for f in rwr_paths2])
    rec_all = np.concatenate([rec1, rec2], axis=0).T

    count1 = np.sum(rec1 > 1.96, axis=0) > rec1.shape[0]/10
    count2 = np.sum(rec2 > 1.96, axis=0) > rec2.shape[0]/10
    rec_sub = rec_all[np.logical_or(count1, count2), :]

    if rec_sub.size == 0:
        return 
    rec_normalized = quantile_normalize(rec_sub, axis=1)
    if any(np.var(rec_normalized, ddof=1, axis=0) < 0):
        print("Variance is negative")
        return 

    loc_info = read_one_rwr_file(rwr_paths1[0], BINSIZE, bids, gids, gap, True)
    loc_info = loc_info.values[np.logical_or(count1, count2), :]*BINSIZE
    out_df = pd.DataFrame(loc_info, columns=["x1", "x2", "y1", "y2"])
    out_df["d"] = out_df["y1"] - out_df["x1"]

    rec_normalized1 = rec_normalized[:, :rec1.shape[0]]
    out_df["mean.A"] = np.mean(rec_normalized1, axis=1)
    rec_normalized2 = rec_normalized[:, rec1.shape[0]:]
    out_df["mean.B"] = np.mean(rec_normalized2, axis=1)

    import warnings
    warnings.filterwarnings("ignore")
    t_stat, p_val = stats.ttest_ind(rec_normalized1, rec_normalized2, axis=1, equal_var=False)
    out_df["Tstat"] = t_stat
    out_df["Ttest.Pvalue"] = p_val
    out_df["fdr"] = multi.multipletests(out_df["Ttest.Pvalue"], method='fdr_bh')[1]
    out_df.insert(0, "chr", chr_id)
    out_df = out_df[np.abs(out_df["mean.A"]-out_df["mean.B"])>1e-10]
    if len(out_df) > 0:
        out_df.to_csv(os.path.join(out_dire, f"test_results_gap{gap}.txt"), sep="\t", index=False)


def call_diff_loops(
        chr_id, dire1, dire2, out_dire, fdr, BINSIZE,
        max_worker, min_gap, max_gap, bid_path, gid_path
):
    if not os.path.exists(out_dire):
        os.mkdir(out_dire)
    bids, gids = get_genome_ids(bid_path, gid_path, BINSIZE, chr_id)
    gaps = np.arange(min_gap, max_gap)
    with ProcessPoolExecutor(max_worker) as executor:
        for result in executor.map(partial(
            single_gap_size, 
            chr_id=chr_id,  BINSIZE=BINSIZE,
            bids=bids, gids=gids, 
            dire1=dire1, dire2=dire2, out_dire=out_dire
        ), gaps):
            pass

    result_fpath = os.path.join(out_dire, "test_results_gap{}.txt")
    out = pd.concat([
        pd.read_csv(result_fpath.format(gap), sep="\t")
        for gap in gaps if os.path.exists(result_fpath.format(gap))
    ])
    out.to_csv(
        os.path.join(out_dire, f"combined_results_{chr_id}.txt"), 
        index=False, sep="\t"
    )
    out[(out["fdr"] < fdr)&(np.abs(out["Tstat"]) > 2)].to_csv(
        os.path.join(out_dire, f"DI_FDR{fdr}_T2_Test_{chr_id}.txt"), 
        index=False, sep="\t"
    )
    for gap in gaps:
        if os.path.exists(result_fpath.format(gap)):
            os.remove(os.path.join(out_dire, f"test_results_gap{gap}.txt"))