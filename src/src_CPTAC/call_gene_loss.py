import numpy as np
import pandas as pd


def call_gene_loss(pairs, samples, cndf):
    print('Calling gene loss...')
    A2_lost = pd.DataFrame(index=pairs.A2.unique(), columns=sorted(samples))  # Will be populate with boolean values indicating gene loss
    cndf['single_copy'] = cndf['CNV'].apply(get_single_copy)
    print('Pivoting and updating...')
    pivot = cndf.pivot_table(index='gene_name', columns='sample_ID', values='single_copy', aggfunc = 'first') # No duplicates, aggfunc specified to speed up
    A2_lost.update(pivot)
    print('Done pivoting and updating!')  # Faster than iterating
    return A2_lost


def get_single_copy(x):  # Hemizygous deletions true, diploid False, rest NaNs. These are GISTIC absolute CN values.
    if x == -1:
        return True
    elif x == 0:
        return False
    else:
        return np.nan