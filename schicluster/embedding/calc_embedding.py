import numpy as np
import pandas as pd
import cooler
from sklearn.decomposition import TruncatedSVD
import joblib
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
import pathlib


def make_idx(n_dim, dist, resolution):
    idx = np.triu_indices(n_dim, k=1)
    # vectorized band filter (Python loop is far too slow/memory-heavy at 25kb)
    idx_filter = (idx[1] - idx[0]) < (dist / resolution + 1)
    idx = (idx[0][idx_filter], idx[1][idx_filter])
    return idx


def make_chrom_matrix_chunk(output_path,
                            cell_urls,
                            row_offset,
                            chrom,
                            idx,
                            scale_factor):
    # multiple workers fill disjoint row ranges of the same pre-created memmap
    # concurrently (rows never overlap, so no locking is needed)
    chrom_matrix = np.lib.format.open_memmap(str(output_path), mode='r+')
    for j, cell_url in enumerate(cell_urls):
        cool = cooler.Cooler(cell_url)
        matrix = cool.matrix(balance=False, sparse=False).fetch(chrom)
        # each row of chrom_matrix is a 1D cell-chrom matrix
        chrom_matrix[row_offset + j, :] = matrix[idx].ravel()
    # cast to float32 first, then scale in float32 to match the original exactly
    chrom_matrix[row_offset:row_offset + len(cell_urls)] *= scale_factor
    chrom_matrix.flush()
    del chrom_matrix
    return


def svd(input_path, dim, output_prefix, save_model=True, norm_sig=True):
    loaded = np.load(input_path)
    # raw chrom matrices are .npy memmaps, the concatenated matrix is a .npz
    chrom_matrix = loaded['arr_0'] if hasattr(loaded, 'files') else loaded
    dim = min(dim, chrom_matrix.shape[0] - 1, chrom_matrix.shape[1] - 1)
    model = TruncatedSVD(n_components=dim, algorithm='arpack')
    decomp = model.fit_transform(chrom_matrix)

    if norm_sig:
        decomp = decomp[:, model.singular_values_ > 0]
        singular_values = model.singular_values_[model.singular_values_ > 0]
        decomp /= singular_values[None, :]

    np.savez(f'{output_prefix}_decomp.npz', decomp)
    if save_model:
        joblib.dump(model, f'{output_prefix}_SVD.lib')
    return


def embedding(cell_table_path,
              output_dir,
              chrom_size_path=None,
              dim=50,
              dist=1000000,
              resolution=100000,
              scale_factor=100000,
              norm_sig=True,
              cpu=1,svd_cpu=1,
              cell_chunk_size=2000,
              save_model=False,
              save_raw=True):
    cell_table = pd.read_csv(cell_table_path,
                             sep='\t',
                             index_col=0,
                             header=None).squeeze(axis=1)
    output_dir = pathlib.Path(output_dir).absolute()
    output_dir.mkdir(parents=True, exist_ok=True)
    raw_dir = output_dir / 'raw'
    raw_dir.mkdir(exist_ok=True)

    first_cool = cooler.Cooler(cell_table.iloc[0])
    chroms = first_cool.chromnames
    chrom_bin_counts = first_cool.bins()[:]['chrom'].value_counts()
    # remove small chroms
    chroms = chrom_bin_counts.index[chrom_bin_counts > 2]
    if chrom_size_path is not None:
        chrom_sizes = pd.read_csv(chrom_size_path, sep='\t', index_col=0, header=None).squeeze(axis=1)
        chroms = chroms.intersection(chrom_sizes.index)

    # prepare raw chromosome 1D matrix
    # Parallelize over cell chunks (not just chromosomes): with hundreds of
    # thousands of cells, per-chromosome tasks would only use ~len(chroms) cores.
    cell_urls = cell_table.tolist()
    n_cells = cell_table.size
    # pre-create one memmap per chrom (and its band index) in the main process,
    # so workers only open it in 'r+' and write disjoint rows
    chrom_idx = {}
    for chrom in chroms:
        nbins = chrom_bin_counts[chrom]
        idx = make_idx(nbins, dist, resolution)
        chrom_idx[chrom] = idx
        mm = np.lib.format.open_memmap(str(raw_dir / f'{chrom}.npy'),
                                       mode='w+', dtype='float32',
                                       shape=(n_cells, idx[0].size))
        del mm

    with ProcessPoolExecutor(cpu) as exe:
        futures = {}
        for chrom in chroms:
            idx = chrom_idx[chrom]
            output_path = raw_dir / f'{chrom}.npy'
            for start in range(0, n_cells, cell_chunk_size):
                end = min(start + cell_chunk_size, n_cells)
                future = exe.submit(make_chrom_matrix_chunk,
                                    output_path,
                                    cell_urls[start:end],
                                    start,
                                    chrom,
                                    idx,
                                    scale_factor)
                futures[future] = (chrom, start)

        for future in as_completed(futures):
            chrom, start = futures[future]
            future.result()
            print(f'{chrom} rows {start}-{start + cell_chunk_size} generated')

    # SVD on each chromosome
    decomp_dir = output_dir / 'decomp'
    decomp_dir.mkdir(exist_ok=True)
    with ProcessPoolExecutor(svd_cpu) as exe:
        futures = {}
        for chrom in chroms:
            chrom_raw_path = raw_dir / f'{chrom}.npy'
            output_prefix = decomp_dir / chrom
            future = exe.submit(svd,
                                input_path=chrom_raw_path,
                                dim=dim,
                                output_prefix=output_prefix,
                                save_model=save_model,
                                norm_sig=norm_sig)
            futures[future] = chrom

        for future in as_completed(futures):
            chrom = futures[future]
            print(f'{chrom} SVD generated')
            future.result()

    # concatenate all chromosome decomp matrix and do another SVD
    total_data = []
    for chrom in chroms:
        decomp_path = decomp_dir / f'{chrom}_decomp.npz'
        data = np.load(str(decomp_path))['arr_0']
        total_data.append(data)
    total_data = np.concatenate(total_data, axis=1)
    total_data_path = decomp_dir / f'total_chrom_decomp_concat.npz'
    np.savez(total_data_path, total_data)

    # final SVD
    output_prefix = decomp_dir / f'total'
    svd(input_path=total_data_path,
        dim=dim,
        output_prefix=output_prefix,
        save_model=save_model,
        norm_sig=norm_sig)

    # clean single chrom
    for chrom in chroms:
        decomp_path = decomp_dir / f'{chrom}_decomp.npz'
        subprocess.run(['rm', '-f', str(decomp_path)])
    if not save_raw:
        subprocess.run(['rm', '-rf', str(raw_dir)])
    return
