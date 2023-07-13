import click
from SnapHiC_D.SnapHiC_D import call_diff_loops

@click.group()
def main():
    pass

@main.command()
@click.option('-i', '--indira', required = True, help = 'input group A directory')
@click.option('-j', '--indirb', required = True, help = 'input group B directory')
@click.option('-o', '--outdir', required = True, help = 'output directory')
@click.option('-c', '--chrom', required = True, help = 'chromosome to process')
@click.option('-b', "--bidpath", required=True, help= 'reference genome filter regions path')
@click.option('-g', "--gidpath", required=True, help= 'referenec genome transcript path')
@click.option('--binsize', type = int, required = False, default = 1e4, help = 'bin size used for binning the reads')
@click.option('-n', '--num_cpus', type = int , required = True, help = 'number of CPUS used for parallel computing')
@click.option('--fdr_threshold', default = 0.1, type = float, required = False, help = 'FDR threshold used for candidate peak detection')
@click.option('--maxi_gap', type = int, default = 100, required = False, help = 'maximum gap between the read pairs')
@click.option('--mini_gap', type = int, default = 2, required = False, help = 'minimum gap between the read pairs')
def diff_loops(
    indira, indirb, outdir, chrom, bidpath, gidpath,
    binsize, num_cpus, fdr_threshold, maxi_gap, mini_gap
):
    call_diff_loops(
        chr_id=chrom, dire1=indira, dire2=indirb, 
        out_dire=outdir, fdr=fdr_threshold, 
        BINSIZE=binsize, max_worker=num_cpus, 
        min_gap=mini_gap, max_gap=maxi_gap, 
        bid_path=bidpath, gid_path=gidpath
    )