:py:mod:`schicluster.impute.snakemake`
======================================

.. py:module:: schicluster.impute.snakemake


Module Contents
---------------

.. py:data:: PACKAGE_DIR

   

.. py:function:: prepare_impute(output_dir, chrom_size_path, output_dist, window_size, step_size, resolution, input_scool=None, cell_table=None, batch_size=100, logscale=False, pad=1, std=1, rp=0.5, tol=0.01, min_cutoff=1e-05, chrom1=1, pos1=2, chrom2=5, pos2=6, cpu_per_job=10)

   prepare snakemake files and directory structure for cell contacts imputation


