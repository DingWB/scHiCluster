import os,sys
import numpy as np
import pandas as pd
import pybedtools
pybedtools.helpers.set_bedtools_path(path=os.path.dirname(sys.executable))
from concurrent.futures import ProcessPoolExecutor, as_completed

def overlap_with_bed(data,bed2,chrom1,pos1,chrom2,pos2):
    data['ID1']=data[chrom1].map(str)+'-'+data[pos1].map(str)
    data['ID2']=data[chrom2].map(str)+'-'+data[pos2].map(str)
    df1=data.loc[:,[chrom1,pos1]].drop_duplicates()
    df1.columns=['chrom','start']
    df1['end']=df1.start + 1
    bed1 = pybedtools.BedTool.from_dataframe(df1)
    intersection = bed1.intersect(bed2, wa=True,wb=True) # loj=True
    df=intersection.to_dataframe(names=df1.columns.tolist()+['chr2','start2','end2','category'])
    if df.shape[0]==0:
        data['category1']='NA'
    else:
        # df.category=df.category.replace({'.':'NA'})
        df['ID']=df.chrom.map(str)+'-'+df.start.map(str)
        D=df.loc[:,['ID','category']].drop_duplicates().set_index('ID').category.to_dict()
        data['category1']=data.ID1.map(D).fillna('NA')

    df1=data.loc[:,[chrom2,pos2]].drop_duplicates()
    df1.columns=['chrom','start']
    df1['end']=df1.start + 1
    bed1 = pybedtools.BedTool.from_dataframe(df1)
    intersection = bed1.intersect(bed2, wa=True,wb=True)
    df=intersection.to_dataframe(names=df1.columns.tolist()+['chr2','start2','end2','category'])
    if df.shape[0]==0:
        data['category2']='NA'
    else:
        df['ID']=df.chrom.map(str)+'-'+df.start.map(str)
        D=df.loc[:,['ID','category']].drop_duplicates().set_index('ID').category.to_dict()
        data['category2']=data.ID2.map(D).fillna('NA')
    data['category']=data.category1.map(str)+'|'+data.category2.map(str)
    # Clean up all temporary files created in the session
    # pybedtools.cleanup()
    return data

def compute_decay(cell_name, contact_path, bins, chrom_sizes=None, resolution=10000, bed_df=None, chrom1=1, chrom2=5, pos1=2, pos2=6):
    # decay
    data = pd.read_csv(contact_path, sep='\t', header=None, index_col=None)
    data = data.loc[(data[chrom1]==data[chrom2]) & data[chrom1].isin(chrom_sizes.index)] # select cis-contact
    if not bed_df is None: # with 4 columns: chrom,start,end and category
        # Intersect the contact dataframe with the given bed regions (such as A, B compartment files)
        try:
            data=overlap_with_bed(data,bed_df,chrom1,pos1,chrom2,pos2)
        except Exception as error_message:
            print(cell_name)
            print(error_message)
        data['dist']=(data[pos2]-data[pos1]).abs()
        data[[pos1, pos2]] = data[[pos1, pos2]] // resolution
        decay_results=[]
        sparsity_results=[]
        for category in data.category.unique():
            data1=data.loc[data.category==category]
            hist=np.histogram(data1.dist.tolist(), bins)[0] # length is 132
            decay_results.append(pd.DataFrame(hist, columns=[f"{cell_name}.{category}"]))
            # sparsity
            data1 = data1.groupby(by=[chrom1, pos1, pos2])[chrom2].count().reset_index() #chrom,pos1,po2s,count
            data1 = data1.loc[data1[pos1]!=data1[pos2], chrom1].value_counts() # total number of cis contacts on each chromosoe
            sparsity_results.append(pd.DataFrame(data1).set_axis([f"{cell_name}.{category}"], axis=1))
        sparsity=pd.concat(sparsity_results,axis=1)
        decay=pd.concat(decay_results,axis=1)
    else:
        hist = np.histogram(np.abs(data[pos2] - data[pos1]), bins)[0] #The number of data points that fall into each bin (the histogram counts).
        decay=pd.DataFrame(hist, columns=[cell_name])
        # sparsity
        data[[pos1, pos2]] = data[[pos1, pos2]] // resolution
        data = data.groupby(by=[chrom1, pos1, pos2])[chrom2].count().reset_index() #chrom,pos1,po2s,count
        data = data.loc[data[pos1]!=data[pos2], chrom1].value_counts() # total number of contacts on each chromosoe
        sparsity=pd.DataFrame(data).set_axis([cell_name], axis=1) # a dataframe of 1 columns (cell_name) if bed_df is None, rows are chromosomes
    return [sparsity,  # sparsity: a dataframe with chroms as index and cell_name as columns names, values are the total number of contacts on each chrom
            decay]  # decay: a dataframe with only one column (column name is cell_name): number of contact fall into each distance bin.

def contact_distance(contact_table=None, chrom_size_path=None, bed_df=None, resolution=10000, output_prefix=None, chrom1=1, chrom2=5, pos1=2, pos2=6, cpu=16,verbose=1):
    outdir=os.path.dirname(output_prefix)
    pybedtools.tempfile.tempdir=outdir
    if not bed_df is None:
        if isinstance(bed_df,str):
            bed_df=pd.read_csv(os.path.expanduser(bed_df),sep='\t',header=None,usecols=[0,1,2],names=['chrom','start','end','category'])
        else:
            assert isinstance(bed_df,pd.DataFrame)
        bed_df = pybedtools.BedTool.from_dataframe(bed_df.drop_duplicates())
    chrom_sizes = pd.read_csv(chrom_size_path, sep='\t', header=None, index_col=0)
    nbins = np.floor(np.log2(chrom_sizes[1].values.max() / 2500) / 0.125) #132
    bins = 2500 * np.exp2(0.125 * np.arange(nbins+1))
    #dist = int(chromsize[1].min() // res + 1)
    if isinstance(contact_table,str):
        contact_table = pd.read_csv(contact_table, sep='\t', index_col=None, header=None)
    else:
        assert isinstance(contact_table,pd.DataFrame) # columns: 0 (cell_name) and 1 (contact_path)
    with ProcessPoolExecutor(cpu) as executor:
        futures = {}
        for cell_name, contact_path in contact_table.values:
            future = executor.submit(
                compute_decay,
                cell_name=cell_name,
                contact_path=contact_path,
                bins=bins,
                chrom_sizes=chrom_sizes,
                resolution=resolution,bed_df=bed_df,
                chrom1=chrom1,
                pos1=pos1,
                chrom2=chrom2,
                pos2=pos2,
            )
            futures[future] = cell_name

        sparsity, decay = [], []
        for future in as_completed(futures):
            cell_name = futures[future]
            xx, yy = future.result()
            sparsity.append(xx)
            decay.append(yy)
            if verbose > 0:
                print(f'{cell_name} finished')
            
    # sparsity = pd.concat(sparsity, axis = 1)[contact_table[0]].T
    # decay = pd.concat(decay, axis = 1)[contact_table[0]].T
    sparsity = pd.concat(sparsity, axis = 1).T
    decay = pd.concat(decay, axis = 1).T
    sparsity.to_hdf(f'{output_prefix}_chromsparsity.hdf5', key='data')
    decay.to_hdf(f'{output_prefix}_decay.hdf5', key='data')
    pybedtools.cleanup(remove_all=True)
