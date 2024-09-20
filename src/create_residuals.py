# cptac_analysis.py

import numpy as np
import pandas as pd
from pandarallel import pandarallel
from tqdm import tqdm
import statsmodels.formula.api as smf
from concurrent.futures import ProcessPoolExecutor
import gc

# Initialize parallel processing
pandarallel.initialize(nb_workers=20, progress_bar=True)

def load_and_preprocess_data():
    """
    Load and preprocess CPTAC data files.
    
    Returns:
    tuple: Processed dataframes (proteomics, transcriptomics, copy number variation)
    """
    # Load data
    cptac_pancan_prot = pd.read_parquet("../data/CPTAC/CPTAC_PanCan_proteomics.parquet")
    cptac_pancan_prot.index.name = 'sample_ID'

    cptac_pancan_trans = pd.read_parquet("../data/CPTAC/CPTAC_PanCan_transcriptomics.parquet")
    cptac_pancan_trans.columns.name = ''
    cptac_pancan_trans.index.name = 'sample_ID'

    cptac_cnv_wide = pd.read_csv('../data/CPTAC/CPTAC_gistic_cn.csv')
    cptac_cnv_wide = cptac_cnv_wide.set_index('gene_name').T
    cptac_cnv_wide.index.name = 'sample_ID'
    cptac_cnv_wide.columns.name = ''
    cptac_cnv_wide = cptac_cnv_wide.parallel_apply(lambda x: pd.to_numeric(x, errors='coerce'))
    cptac_cnv_wide = np.floor(cptac_cnv_wide).astype('Int64')

    # Get common genes and samples
    common_genes = list(set(cptac_pancan_prot.columns) & set(cptac_pancan_trans.columns) & set(cptac_cnv_wide.columns))
    common_samples = list(set(cptac_pancan_prot.index) & set(cptac_pancan_trans.index) & set(cptac_cnv_wide.index))

    # Filter and sort data
    cptac_cnv_wide = cptac_cnv_wide.astype(float)

    cptac_pancan_prot = cptac_pancan_prot.loc[common_samples, common_genes].sort_index().sort_index(axis=1)
    cptac_pancan_trans = cptac_pancan_trans.loc[common_samples, common_genes].sort_index().sort_index(axis=1)
    cptac_cnv_wide = cptac_cnv_wide.loc[common_samples, common_genes].sort_index().sort_index(axis=1)

    return cptac_pancan_prot, cptac_pancan_trans, cptac_cnv_wide

def melt_dataframe(df, value_name):
    """
    Melt dataframe to long format.
    
    Args:
    df (pd.DataFrame): Input dataframe
    value_name (str): Name for the value column
    
    Returns:
    pd.DataFrame: Melted dataframe
    """
    melted_df = df.reset_index().melt(id_vars='sample_ID', var_name='gene_name', value_name=value_name)
    melted_df['gene_in_sample'] = melted_df['gene_name'] + '_in_' + melted_df['sample_ID']
    return melted_df.drop(columns=['sample_ID', 'gene_name'])

def merge_dataframes(prot_long, trans_long, cnv_long):
    """
    Merge proteomics, transcriptomics, and copy number variation dataframes.
    
    Args:
    prot_long (pd.DataFrame): Proteomics dataframe
    trans_long (pd.DataFrame): Transcriptomics dataframe
    cnv_long (pd.DataFrame): Copy number variation dataframe
    
    Returns:
    pd.DataFrame: Merged dataframe
    """
    merged_df = prot_long.merge(trans_long, on='gene_in_sample').merge(cnv_long, on='gene_in_sample')
    merged_df['sample_ID'] = merged_df['gene_in_sample'].apply(lambda x: x.split('_in_')[1])

    # Load and merge sample info
    cptac_sample_info = pd.read_csv('../data/CPTAC/CPTAC_sampleinfo.csv').rename(columns={'Unnamed: 0': 'sample_ID'})
    cptac_sample_info['sample_ID'] = cptac_sample_info['sample_ID'].str.lstrip('X')

    merged_df = merged_df.merge(cptac_sample_info[['sample_ID', 'Study']], on='sample_ID', how='inner')
    merged_df = merged_df[merged_df.Study != 'OvCa 2020']
    return merged_df.drop(columns='sample_ID')

def process_gene(gene_data):
    gene, gene_df, dependent_var, independent_vars = gene_data
    formula = f'{dependent_var} ~ ' + ' + '.join(independent_vars)
    model = smf.ols(formula, data=gene_df).fit()
    residuals = model.resid
    return gene, gene_df, residuals

def regress_out_variables_optimized(df, dependent_var, independent_vars, n_workers=20):
    print("Preparing data for regression...")
    df[['Gene', 'Sample']] = df['gene_in_sample'].str.split('_in_', expand=True)
    
    print("Grouping data by gene...")
    gene_groups = df.groupby('Gene')
    
    print("Preparing gene data for parallel processing...")
    gene_data = []
    for gene, group in gene_groups:
        gene_df = group[['Sample', dependent_var] + independent_vars].dropna()
        if len(gene_df) > 1:
            gene_data.append((gene, gene_df, dependent_var, independent_vars))
    
    total_genes = len(gene_data)
    print(f"Prepared {total_genes} genes for processing")
    
    residuals_dict = {}
    
    print(f"Starting parallel processing of {total_genes} genes...")
    
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        results = list(tqdm(executor.map(process_gene, gene_data), total=total_genes, desc="Genes processed"))
        
        for gene, gene_df, residuals in results:
            residuals_dict.update({f"{gene}_in_{sample}": residual for sample, residual in zip(gene_df['Sample'], residuals)})
    
    print(f"Finished processing {total_genes} genes.")
    
    print("Creating final dataframe...")
    residuals_df = pd.DataFrame(list(residuals_dict.items()), columns=['gene_in_sample', 'Residuals'])    
    print("Regression complete.")
    return residuals_df

def main():
    # Load and preprocess data
    print("Loading and preprocessing data...")
    cptac_pancan_prot, cptac_pancan_trans, cptac_cnv_wide = load_and_preprocess_data()

    # Melt dataframes to long format
    print("Melting dataframes...")
    cptac_pancan_prot_long = melt_dataframe(cptac_pancan_prot, 'prot_abundance')
    cptac_pancan_trans_long = melt_dataframe(cptac_pancan_trans, 'trans_abundance')
    cptac_cnv_long = melt_dataframe(cptac_cnv_wide, 'copy_number')

    # Merge dataframes
    print("Merging dataframes...")
    pred_df = merge_dataframes(cptac_pancan_prot_long, cptac_pancan_trans_long, cptac_cnv_long)

    print('Protein minus study + CN running')
    prot_minus_study_CN = regress_out_variables_optimized(pred_df, 'prot_abundance', ['Study', 'copy_number'])
    prot_minus_study_CN.to_csv('../../modular3_newresids/prot_study_CN_resid.csv', index=False)
    gc.collect()

    print('Transcript minus study + CN running')
    trans_minus_study_CN = regress_out_variables_optimized(pred_df, 'trans_abundance', ['Study', 'copy_number'])
    trans_minus_study_CN.to_csv('../../modular3_newresids/trans_study_CN_resid.csv', index=False)
    gc.collect()

    # Regress out study effect from protein abundance
    print('Protein minus study running')
    prot_minus_study = regress_out_variables_optimized(pred_df, 'prot_abundance', ['Study'])
    prot_minus_study.to_csv('../../modular3_newresids/prot_study_resid.csv', index=False)
    gc.collect()
   
    # Regress out study effect from transcript abundance
    print('Transcript minus study running')
    trans_minus_study = regress_out_variables_optimized(pred_df, 'trans_abundance', ['Study'])
    trans_minus_study.to_csv('../../modular3_newresids/trans_study_resid.csv', index=False)
    gc.collect()

    # Regress out study-corrected transcriptomics from study-corrected proteomics
    print('Protein residual calculation running')
    prot_residual = regress_out_variables_optimized(
    prot_minus_study.rename(columns={'Residuals': 'Residuals_prot_study'}).merge(
    trans_minus_study[['gene_in_sample', 'Residuals']].rename(columns={'Residuals': 'Residuals_trans_study'}),
            on='gene_in_sample'),
        'Residuals_prot_study',  ['Residuals_trans_study'])
    prot_residual.to_csv('../../modular3_newresids/CPTAC_prot_residuals.csv')
    gc.collect()

if __name__ == "__main__":
    main()