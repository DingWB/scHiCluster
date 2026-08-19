import cooler
import numpy as np
import pandas as pd
from scipy.sparse import triu, csr_matrix
from concurrent.futures import ProcessPoolExecutor, as_completed


def _oe_decay(D, max_off):
    """Mean contact at each diagonal offset 1..max_off (distance expected).

    Used for the observed/expected gene score. D is an upper-triangular sparse
    matrix; each diagonal's mean uses its full length (n - offset) so unimputed
    zeros count, matching the standard distance-decay definition.
    """
    n = D.shape[0]
    coo = D.tocoo()
    offset = coo.col - coo.row
    keep = (offset >= 1) & (offset <= max_off)
    sums = np.bincount(offset[keep], weights=coo.data[keep], minlength=max_off + 1)
    counts = np.maximum(n - np.arange(max_off + 1), 1)
    return sums / counts


def _window_score(data, offset, decay, oe_normalize, per_gene_mean, eps=1e-5):
    """Reduce one gene window's contacts to a scalar score.

    oe_normalize divides each contact by the distance-decay expectation at its
    offset (observed/expected, removing distance bias); per_gene_mean divides by
    the number of contacts in the window (removing the gene-length effect).
    """
    if data.size == 0:
        return 0.0
    vals = data / (decay[offset] + eps) if oe_normalize else data
    s = float(vals.sum())
    if per_gene_mean:
        s = s / data.size
    return s


def gene_score_raw(cell_path, chrom_sizes, gene_meta, resolution, chrom1, pos1,
                   chrom2, pos2, per_gene_mean=False):
    data = pd.read_csv(cell_path, sep='\t', index_col=None, header=None, comment='#')
    data = data.loc[(data[chrom1]==data[chrom2]) & data[chrom1].isin(chrom_sizes.index)]
    result = []
    for chrom in chrom_sizes.index:
        n_bins = int(chrom_sizes.loc[chrom, 1] // resolution) + 1
        chrfilter = (data[chrom1]==chrom)
        if chrfilter.sum()==0:
            D = csr_matrix((n_bins, n_bins))
        else:
            D = data.loc[chrfilter].copy()
            D[[pos1, pos2]] = (D[[pos1, pos2]] - 1) // resolution
            D = D.groupby(by=[pos1, pos2])[chrom1].count().reset_index()
            D = csr_matrix((D[chrom1].astype(np.int32), (D[pos1], D[pos2])), (n_bins, n_bins))
        gene = gene_meta.loc[gene_meta[0]==chrom, [1,2]].values
        for xx,yy in gene:
            sub = D[xx:(yy+1), xx:(yy+1)]
            if per_gene_mean:
                sub = sub.tocoo()
                result.append(float(sub.data.sum()) / sub.nnz if sub.nnz else 0.0)
            else:
                result.append(sub.sum())
    return result

def gene_score_impute(cell_path, chrom_sizes, gene_meta,
                      oe_normalize=False, per_gene_mean=False):
    cool = cooler.Cooler(cell_path)
    result = []
    for chrom in chrom_sizes.index:
        D = triu(cool.matrix(balance=False, sparse=True).fetch(chrom), k=1).tocsr()
        gene = gene_meta.loc[gene_meta[0]==chrom, [1,2]].values
        decay = None
        if oe_normalize and gene.shape[0] > 0 and D.shape[0] > 1:
            # distance-decay expected per offset, capped at the widest gene window
            max_off = min(int(gene[:, 1].max() - gene[:, 0].min()) + 2, D.shape[0] - 1)
            decay = _oe_decay(D, max_off)
        for xx,yy in gene:
            sub = D[(xx-1):(yy+1), xx:(yy+2)].tocoo()
            if not oe_normalize and not per_gene_mean:
                result.append(float(sub.data.sum()))
                continue
            # window offset = abs_col - abs_row = (local col - local row) + 1
            offset = sub.col - sub.row + 1
            result.append(_window_score(sub.data, offset, decay,
                                        oe_normalize, per_gene_mean))
    return result

def gene_score(cell_table_path, gene_meta_path, resolution, output_hdf_path, chrom_size_path, 
               slop=0, cpu=10, mode='impute', chrom1=1, pos1=2, chrom2=5, pos2=6,
               oe_normalize=False, per_gene_mean=False, min_gene_bins=None):
    chrom_sizes = pd.read_csv(chrom_size_path, sep='\t', header=None, index_col=0)
    gene_meta = pd.read_csv(gene_meta_path, sep='\t', header=None, index_col=3)
    gene_meta = gene_meta[gene_meta[0].isin(chrom_sizes.index)]
    gene_meta[1] = (gene_meta[1] - slop) // resolution
    gene_meta[2] = (gene_meta[2] + slop) // resolution
    # (C) widen very short genes to at least min_gene_bins bins (symmetric) to
    # stabilize their score; clip start to >=1 so the (xx-1) impute window is valid.
    if min_gene_bins and min_gene_bins > 1:
        span = gene_meta[2] - gene_meta[1]
        need = ((min_gene_bins - 1) - span).clip(lower=0)
        left = need // 2
        gene_meta[1] = (gene_meta[1] - left).clip(lower=1)
        gene_meta[2] = gene_meta[2] + (need - left)
    if mode == 'raw' and oe_normalize:
        print("WARNING: oe_normalize is only supported for mode='impute'; ignored for raw.")
    # the workers append gene results grouped by chrom in chrom_sizes.index order,
    # not in gene_meta file order, so the output columns must follow that same
    # order -- otherwise gene ids get misaligned with their values.
    gene_order = [gid for chrom in chrom_sizes.index
                  for gid in gene_meta.index[gene_meta[0] == chrom]]
    cell_table = pd.read_csv(cell_table_path, sep='\t', header=None, index_col=None).values
    with ProcessPoolExecutor(cpu) as exe:
        future_dict = {}
        for cell, cell_path in cell_table:
            if mode=='impute':
                future = exe.submit(gene_score_impute,
                                    cell_path=cell_path,
                                    chrom_sizes=chrom_sizes,
                                    gene_meta=gene_meta,
                                    oe_normalize=oe_normalize,
                                    per_gene_mean=per_gene_mean)
            elif mode=='raw':
                future = exe.submit(gene_score_raw,
                                    cell_path=cell_path,
                                    chrom_sizes=chrom_sizes,
                                    gene_meta=gene_meta, 
                                    resolution=resolution, 
                                    chrom1=chrom1, 
                                    pos1=pos1, 
                                    chrom2=chrom2, 
                                    pos2=pos2,
                                    per_gene_mean=per_gene_mean)
            else:
                print("ERROR : Mode must be one of impute/raw/diag")
                return 0
            future_dict[future] = cell

        result, cell_list = [], []
        for future in as_completed(future_dict):
            cell = future_dict[future]
            print(f'{cell} finished.')
            result.append(future.result())
            cell_list.append(cell)

    result = pd.DataFrame(result, index=cell_list, columns=gene_order)
    # restore the original gene_meta order for reproducibility (labels are now
    # correctly aligned, so this is a safe relabeling-free reorder)
    result = result.reindex(columns=gene_meta.index)
    result.to_hdf(output_hdf_path, key='data', complevel=9)

