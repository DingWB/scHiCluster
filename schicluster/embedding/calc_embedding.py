import numpy as np
import pandas as pd
import cooler
from cooler.util import parse_cooler_uri
import h5py
from sklearn.utils.extmath import randomized_svd
import joblib
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
import pathlib


def make_idx(n_dim, dist, resolution):
    # kept for reference / backwards compatibility. The fill path below uses
    # band_row_offsets() instead, which is O(n_dim) rather than O(n_dim**2):
    # np.triu_indices(9959) alone allocates ~800 MB before the band filter.
    idx = np.triu_indices(n_dim, k=1)
    idx_filter = (idx[1] - idx[0]) < (dist / resolution + 1)
    idx = (idx[0][idx_filter], idx[1][idx_filter])
    return idx


def band_max_offset(dist, resolution):
    # largest diagonal offset make_idx() keeps, i.e. the largest integer
    # strictly below dist / resolution + 1
    return int(np.ceil(dist / resolution + 1)) - 1


def band_row_offsets(n_dim, max_off):
    """P[i] = position in the flat band vector where row i starts.

    Row i holds columns i+1 .. min(i + max_off, n_dim - 1), which is exactly
    the set make_idx() keeps, in the same row-major order -- so the column
    layout of the output matrices is unchanged. P[-1] is the band length.
    """
    rows = np.arange(n_dim)
    lens = np.minimum(rows + max_off, n_dim - 1) - rows
    P = np.zeros(n_dim + 1, dtype=np.int64)
    np.cumsum(lens, out=P[1:])
    return P


def band_offset_index(row_offsets):
    """Diagonal offset (1..max_off) of each position in the flat band vector.

    Position p in row i (see band_row_offsets) holds M[i, i+offset], so
    offset = p - row_offsets[i] + 1. Used to group band entries by genomic
    distance for the observed/expected normalization.
    """
    total = int(row_offsets[-1])
    idx = np.arange(total)
    row_id = np.searchsorted(row_offsets, idx, side='right') - 1
    return (idx - row_offsets[row_id] + 1).astype(np.int32)


def _normalize_band(band, off_index, off_counts, max_off, scale,
                    distance_normalize, log_transform, eps=np.float32(1e-5)):
    """Optionally distance-normalize (O/E) and log-transform one band vector.

    distance_normalize divides each entry by the mean contact at its diagonal
    offset (observed/expected), which removes the distance decay that otherwise
    lets the first few diagonals dominate the SVD and cancels most of the
    per-cell depth. log_transform then applies log1p. With both disabled it
    falls back to the original ``band * scale`` behaviour.
    """
    if distance_normalize:
        sums = np.bincount(off_index, weights=band.astype(np.float64),
                           minlength=max_off + 1)
        mean_per_off = (sums / off_counts).astype(np.float32)
        band = band / (mean_per_off[off_index] + eps)
    else:
        band = band * scale
    if log_transform:
        band = np.log1p(band)
    return band.astype(np.float32, copy=False)


def _cell_chrom_band(grp, bin1_offset, first_bin, last_bin, row_offsets, max_off):
    """Read one cell/chromosome band vector straight from the pixel table.

    The old path called cool.matrix(sparse=False).fetch(chrom), which densifies
    the whole chromosome (793 MB for chr1 at 25kb) just to keep the ~40
    diagonals within `dist`. Benchmarking showed fetch() splits about evenly
    between gzip decompression and cooler's matrix assembly; sparse=True only
    removes the last step of the latter and is not faster overall. Going
    straight from the pixel slices to the band vector removes the assembly
    entirely: 2.44x faster per cell and ~25x lower peak memory, bit-identical
    output.

    Requires storage-mode 'symmetric-upper' (checked in embedding()), so
    bin2_id >= bin1_id and each band entry appears exactly once.
    """
    lo = int(bin1_offset[first_bin])
    hi = int(bin1_offset[last_bin])
    band = np.zeros(int(row_offsets[-1]), dtype=np.float32)
    if hi <= lo:
        return band

    bin1 = grp['pixels/bin1_id'][lo:hi]
    bin2 = grp['pixels/bin2_id'][lo:hi]
    offset = bin2 - bin1
    # bin2 < last_bin also drops any trans pixel that would otherwise land in
    # this chromosome's band at the chromosome boundary
    keep = (offset >= 1) & (offset <= max_off) & (bin2 < last_bin)
    if not keep.any():
        return band

    count = grp['pixels/count'][lo:hi]
    # imputed coolers can store NaN/inf counts (row/col normalization over a
    # zero-coverage bin); copied verbatim they poison the whole offset's O/E mean
    count = np.nan_to_num(count, nan=0.0, posinf=0.0, neginf=0.0)
    band[row_offsets[bin1[keep] - first_bin] + (offset[keep] - 1)] = count[keep]
    return band


def make_chrom_matrix_chunk(output_paths,
                            cell_urls,
                            row_offset,
                            chroms,
                            chrom_band,
                            scale_factor,
                            max_off,
                            distance_normalize=True,
                            log_transform=True):
    # scale_factor is only applied on the legacy path (distance_normalize=False);
    # when distance_normalize is on the band is O/E-normalized and scale is unused.
    # each worker handles one disjoint range of cells for ALL chroms at once,
    # so every cool file is opened only once (instead of once per chrom).
    # workers fill disjoint row ranges of the pre-created memmaps concurrently
    # (rows never overlap, so no locking is needed).
    chrom_matrices = {chrom: np.lib.format.open_memmap(str(output_paths[chrom]), mode='r+')
                      for chrom in chroms}
    scale = np.float32(scale_factor)
    for j, cell_url in enumerate(cell_urls):
        store_path, group_path = parse_cooler_uri(cell_url)
        with h5py.File(store_path, 'r') as store:
            grp = store[group_path]
            bin1_offset = grp['indexes/bin1_offset'][:]
            for chrom in chroms:
                first_bin, last_bin, row_offsets, off_index, off_counts = chrom_band[chrom]
                band = _cell_chrom_band(grp, bin1_offset, first_bin, last_bin,
                                        row_offsets, max_off)
                # distance-normalize (O/E) + log1p in memory; the old code only
                # did band *= scale then a second pass through the memmap
                band = _normalize_band(band, off_index, off_counts, max_off,
                                       scale, distance_normalize, log_transform)
                chrom_matrices[chrom][row_offset + j, :] = band
    for chrom in chroms:
        chrom_matrices[chrom].flush()
    del chrom_matrices
    return


def svd(input_path, dim, output_prefix, save_model=True, norm_sig=True,
        n_iter=4, n_oversamples=150, fit_n_cells=None):
    # mmap_mode keeps the raw chrom matrix on disk (chr1 is 364 GB at 25kb for
    # 228k cells); it is silently ignored for the .npz concat matrix, which is
    # small enough to load whole.
    loaded = np.load(input_path, mmap_mode='r')
    # raw chrom matrices are .npy memmaps, the concatenated matrix is a .npz
    chrom_matrix = loaded['arr_0'] if hasattr(loaded, 'files') else loaded
    n_rows, n_features = chrom_matrix.shape
    dim = min(dim, n_rows - 1, n_features - 1)
    # randomized_svd rather than TruncatedSVD: TruncatedSVD always computes
    # explained_variance_ratio_ via np.var(X, axis=0), which materialises a
    # full-size temporary and defeats the memmap entirely. randomized_svd
    # streams the memmap and only ever allocates n_samples x (dim +
    # n_oversamples). arpack is not an option at this scale -- it needs
    # hundreds of passes over a 364 GB matrix.
    #
    # These band matrices have a very flat tail (on chr1, S[48]/S[49] = 1.002),
    # so power iterations have almost no spectral gap to amplify and
    # oversampling is what buys accuracy. n_oversamples is nearly free: it
    # widens the sketch but does not add passes over the data.

    if fit_n_cells is not None and fit_n_cells < n_rows:
        # The SVD "model" (components + singular values) lives in feature space
        # and is independent of the cell count, so fit it on a random cell
        # subset and then project every cell onto it: X @ components.T == U * S.
        # This is one streaming pass instead of the (2*n_iter+2) passes a full
        # randomized_svd makes. Subset peak memory is fit_n_cells*n_features*4B.
        rng = np.random.default_rng(0)
        sub_idx = np.sort(rng.choice(n_rows, size=fit_n_cells, replace=False))
        fit_matrix = np.asarray(chrom_matrix[sub_idx])
        n_os = max(10, min(n_oversamples, min(fit_matrix.shape) - dim))
        _, singular_values, components = randomized_svd(
            fit_matrix, n_components=dim, n_oversamples=n_os,
            n_iter=n_iter, random_state=0)
        del fit_matrix
        # project all cells block-wise so the full matrix is never in RAM
        Vt = np.ascontiguousarray(components.T, dtype=np.float32)
        decomp = np.empty((n_rows, dim), dtype=np.float32)
        block = max(1, int(5e8 // (n_features * 4)))
        for start in range(0, n_rows, block):
            end = min(start + block, n_rows)
            decomp[start:end] = np.asarray(chrom_matrix[start:end]) @ Vt
    else:
        n_os = max(10, min(n_oversamples, min(chrom_matrix.shape) - dim))
        U, singular_values, components = randomized_svd(
            chrom_matrix, n_components=dim, n_oversamples=n_os,
            n_iter=n_iter, random_state=0)
        decomp = U * singular_values

    if norm_sig:
        keep = singular_values > 0
        decomp = decomp[:, keep] / singular_values[keep][None, :]

    np.savez(f'{output_prefix}_decomp.npz', decomp)
    if save_model:
        # randomized_svd returns arrays, not a fitted estimator, so this is a
        # dict rather than the sklearn TruncatedSVD object the old code dumped
        joblib.dump({'components': components,
                     'singular_values': singular_values,
                     'n_components': dim},
                    f'{output_prefix}_SVD.lib')
    return


def embedding(cell_table_path,
              output_dir,
              chrom_size_path=None,
              dim=50,
              dist=1000000,
              resolution=100000,
              scale_factor=100000,
              norm_sig=True,
              cpu=24, svd_cpu=4,
              cell_chunk_size=500,
              save_model=False,
              save_raw=True,
              svd_n_iter=4,
              svd_n_oversamples=150,
              svd_downsample=100000,
              distance_normalize=False,
              log_transform=False):
    cell_table = pd.read_csv(cell_table_path,
                             sep='\t',
                             index_col=0,
                             header=None).squeeze(axis=1)
    output_dir = pathlib.Path(output_dir).absolute()
    output_dir.mkdir(parents=True, exist_ok=True)
    raw_dir = output_dir / 'raw'
    raw_dir.mkdir(exist_ok=True)

    first_cool = cooler.Cooler(cell_table.iloc[0])
    if first_cool.storage_mode != 'symmetric-upper':
        raise ValueError(
            f'the pixel-table fill path requires storage-mode '
            f'"symmetric-upper", got "{first_cool.storage_mode}"')
    chroms = first_cool.chromnames
    chrom_bin_counts = first_cool.bins()[:]['chrom'].value_counts()
    # remove small chroms
    chroms = chrom_bin_counts.index[chrom_bin_counts > 2]
    if chrom_size_path is not None:
        chrom_sizes = pd.read_csv(chrom_size_path, sep='\t', index_col=0, header=None).squeeze(axis=1)
        chroms = chroms.intersection(chrom_sizes.index)

    # resolve the raw matrix path for each chrom, accepting either the .npy
    # memmap this function writes or a pre-existing .npz (svd() reads both)
    def resolve_raw_path(chrom):
        npy_path = raw_dir / f'{chrom}.npy'
        npz_path = raw_dir / f'{chrom}.npz'
        if npy_path.exists():
            return npy_path
        if npz_path.exists():
            return npz_path
        return None

    # skip the (expensive) raw generation when every chrom already has its
    # per-chromosome raw file on disk, and go straight to SVD
    raw_paths = {chrom: resolve_raw_path(chrom) for chrom in chroms}
    raw_exists = all(path is not None for path in raw_paths.values())

    if raw_exists:
        print(f'raw matrices already exist in {raw_dir}, '
              f'skipping raw generation')
    else:
        # prepare raw chromosome 1D matrix
        # Parallelize over cell chunks (not just chromosomes): with hundreds of
        # thousands of cells, per-chromosome tasks would only use ~len(chroms) cores.
        cell_urls = cell_table.tolist()
        n_cells = cell_table.size
        max_off = band_max_offset(dist, resolution)
        # pre-create one memmap per chrom (and its band row offsets) in the main
        # process, so workers only open it in 'r+' and write disjoint rows.
        # chrom_band holds ~1 MB total; the old chrom_idx held the full triu index
        # pair (79 MB at 25kb), which was pickled to every worker task.
        chrom_band = {}
        output_paths = {}
        for chrom in chroms:
            first_bin, last_bin = first_cool.extent(chrom)
            nbins = last_bin - first_bin
            row_offsets = band_row_offsets(nbins, max_off)
            off_index = band_offset_index(row_offsets)
            off_counts = np.maximum(
                np.bincount(off_index, minlength=max_off + 1), 1).astype(np.float64)
            chrom_band[chrom] = (first_bin, last_bin, row_offsets,
                                 off_index, off_counts)
            output_paths[chrom] = raw_dir / f'{chrom}.npy'
            mm = np.lib.format.open_memmap(str(output_paths[chrom]),
                                           mode='w+', dtype='float32',
                                           shape=(n_cells, int(row_offsets[-1])))
            del mm
        # guard the O(n) band layout against the original O(n^2) one, on the
        # smallest chrom so make_idx stays cheap
        check_chrom = min(chroms, key=lambda c: chrom_bin_counts[c])
        assert int(chrom_band[check_chrom][2][-1]) == \
            make_idx(int(chrom_bin_counts[check_chrom]), dist, resolution)[0].size

        with ProcessPoolExecutor(cpu) as exe:
            futures = {}
            for start in range(0, n_cells, cell_chunk_size):
                end = min(start + cell_chunk_size, n_cells)
                future = exe.submit(make_chrom_matrix_chunk,
                                    output_paths,
                                    cell_urls[start:end],
                                    start,
                                    list(chroms),
                                    chrom_band,
                                    scale_factor,
                                    max_off,
                                    distance_normalize,
                                    log_transform)
                futures[future] = start

            for future in as_completed(futures):
                start = futures[future]
                future.result()
                print(f'cell rows {start}-{start + cell_chunk_size} generated')

        raw_paths = {chrom: output_paths[chrom] for chrom in chroms}

    # SVD on each chromosome
    # downsample only the per-chrom SVD (the expensive step); <=0 disables it
    fit_n_cells = svd_downsample if svd_downsample and svd_downsample > 0 else None
    decomp_dir = output_dir / 'decomp'
    decomp_dir.mkdir(exist_ok=True)
    with ProcessPoolExecutor(svd_cpu) as exe:
        futures = {}
        for chrom in chroms:
            # skip chroms whose per-chromosome decomp already exists
            if (decomp_dir / f'{chrom}_decomp.npz').exists():
                print(f'{chrom} decomp already exists, skipping SVD')
                continue
            chrom_raw_path = raw_paths[chrom]
            output_prefix = decomp_dir / chrom
            future = exe.submit(svd,
                                input_path=chrom_raw_path,
                                dim=dim,
                                output_prefix=output_prefix,
                                save_model=save_model,
                                norm_sig=norm_sig,
                                n_iter=svd_n_iter,
                                n_oversamples=svd_n_oversamples,
                                fit_n_cells=fit_n_cells)
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

    # final SVD on the concat matrix (small n_features) is left at full
    # resolution so the embedding actually used downstream is exact
    output_prefix = decomp_dir / f'total'
    svd(input_path=total_data_path,
        dim=dim,
        output_prefix=output_prefix,
        save_model=save_model,
        norm_sig=norm_sig,
        n_iter=svd_n_iter,
        n_oversamples=svd_n_oversamples)

    # clean single chrom
    for chrom in chroms:
        decomp_path = decomp_dir / f'{chrom}_decomp.npz'
        subprocess.run(['rm', '-f', str(decomp_path)])
    if not save_raw:
        subprocess.run(['rm', '-rf', str(raw_dir)])
    return
