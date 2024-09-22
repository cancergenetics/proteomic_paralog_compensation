import numpy as np
import pandas as pd
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests
from data_preprocessing import load_and_process_data
from generate_overlaps import generate_overlaps
from copy import deepcopy
import csv


def map_uniprot_to_symbol(uniprot, mapping):
    ids = uniprot.split(";")
    res = ""
    for uniprot_id in ids:
        if uniprot_id in mapping:
            if res:
                res = res + ";" + mapping[uniprot_id]
            else:
                res = mapping[uniprot_id]
    if not res:
        res = f"UNKNOWN_{uniprot}"
    return res


def map_uniprot_to_entrez(uniprot, mapping):
    ids = uniprot.split(";")
    res = ""
    for uniprot_id in ids:
        if uniprot_id in mapping:
            if res:
                res = res + ";" + mapping[uniprot_id]
            else:
                res = mapping[uniprot_id]
    if not res:
        res = f"UNKNOWN_{uniprot}"
    return res


def process_proteomics_data(
    filename,
    symbol_map,
    entrez_map,
    dropped_samples,
    percent_present=0.8,
    fill=True,
    drop_missing=True,
    log_transform=True,
    adjust_means=True
):
    data = pd.read_csv(filename, sep="\t", skiprows=[1])
    data = data[data['Reverse'] != '+']
    data = data[data['Potential contaminant'] != '+']
    description_columns = list(data.columns[:15])
    lfq_columns = [x for x in data.columns if 'LFQ' in x and x not in dropped_samples]
    short_data = data[description_columns + lfq_columns]
    short_data[lfq_columns] = short_data[lfq_columns].replace(to_replace=0, value=np.NAN)
    short_data[lfq_columns] = short_data[lfq_columns].replace(to_replace=1, value=np.NAN)
    if log_transform:
        short_data[lfq_columns] = short_data[lfq_columns].apply(np.log2)
    if adjust_means:
        full_dataset = short_data[lfq_columns].dropna(how='any')
        adjustment = full_dataset.mean() - full_dataset.mean().mean()
        short_data[lfq_columns] = short_data[lfq_columns].sub(adjustment)
    if drop_missing:
        short_data = short_data.dropna(
            axis=0,
            thresh=len(lfq_columns) * percent_present,
            subset=list(lfq_columns)
        )
    if fill:
        short_data[lfq_columns] = short_data[lfq_columns].fillna(short_data[lfq_columns].min())
    short_data["GeneNames"] = short_data["Majority protein IDs"].apply(
        map_uniprot_to_symbol, args=(symbol_map,)
    )
    short_data["EntrezIDs"] = short_data["Majority protein IDs"].apply(
        map_uniprot_to_entrez, args=(entrez_map,)
    )
    return short_data


def run_ttest_for_clone(gene, df):
    gene_old = deepcopy(gene)
    gene = gene.split('_')[0]
    genelist = [x for x in df.index.to_list() if ((f';{gene}' in x) or (f'{gene};' in x) or x == gene)]
    if len(genelist) > 0:
        gene_index = genelist[0]
        KO_quants = df[(df.index == gene_index)][[x for x in df.columns.to_list() if gene_old in x]].iloc[0].values.tolist()
        WT_quants = df[(df.index == gene_index)][[x for x in df.columns.to_list() if 'WT_' in x]].iloc[0].values.tolist()
        return stats.ttest_ind(KO_quants, WT_quants, equal_var=True)
    else:
        return [np.nan, np.nan]


def calculate_logFC(data, key):
    if key in list(data.keys()):
        ko_data, wt_data = data[key]
        ko_mean = np.mean(ko_data)
        wt_mean = np.mean(wt_data)
        if ko_mean <= 0 or wt_mean <= 0:
            return np.nan
        logFC = ko_mean - wt_mean
        return logFC
    else:
        return np.nan


def drop_lfq(inputstr):
    if "_KO_" in inputstr:
        split_str = inputstr.split("_KO_")
        split_str2 = split_str[0].split(" ")[-1]
        return f'{split_str2}_{split_str[1]}'
    else:
        return inputstr.split(" ")[-1]


def group_A2_clones(df):
    df['clone'] = df['A2'].apply(lambda x: x.split('_')[1] if len(x.split('_')) > 1 else None)
    df['real_A2'] = df['A2'].apply(lambda x: x.split('_')[0])
    grouped = df.groupby('real_A2')

    def select_clone(group):
        if len(group) == 1:
            return group.iloc[0]  # If only one row, return it as is
        sig_clones = group[group['sig'] == True]
        if len(sig_clones) == 0:
            return group.iloc[0]  # If no significant clones, return the first row
        elif len(sig_clones) == 1:
            return sig_clones.iloc[0]  # If only one significant clone, return it
        else:
            # If multiple significant clones, return the one with highest absolute logFC
            return sig_clones.loc[sig_clones['logFC'].abs().idxmax()]

    result = grouped.apply(select_clone)
    result = result.reset_index(drop=True)
    return result


def run_ttest_for_paralog(genepair, A2_df_input, A2_df_input_agg, df):
    paralog = genepair.split('_')[0]
    gene = genepair.split('_')[1]
    if len(A2_df_input[A2_df_input.A2 == 'gene'] > 1):
        clone = A2_df_input_agg.set_index('A2').loc[gene, 'clone']
        gene = f'{gene}_{clone}'
    genelist = [x for x in df.index.to_list() if ((f';{paralog}' in x) or (f'{paralog};' in x) or (x == paralog))]
    if len(genelist) > 0:
        gene_index = genelist[0]
        KO_quants = df[(df.index == gene_index)][[x for x in df.columns.to_list() if gene in x]].iloc[0].values.tolist()
        WT_quants = df[(df.index == gene_index)][[x for x in df.columns.to_list() if 'WT_' in x]].iloc[0].values.tolist()
        return stats.ttest_ind(KO_quants, WT_quants, equal_var=True)
    else:
        return [np.nan, np.nan]


def main():

    print('Processing HAP1 proteomic data and running tests...')

    filename_input = '../data/HAP1/proteinGroups_rerunaug.txt'
    uniprot_to_symbol = {}
    uniprot_to_entrez = {}
    with open("../data/general/geneidmap_sep24.txt", 'r') as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            if r['Approved symbol'] and r['UniProt ID(supplied by UniProt)']:
                parts = r['UniProt ID(supplied by UniProt)'].split(';')
                for p in parts:
                    uniprot_to_symbol[p] = r['Approved symbol']
            if r['NCBI Gene ID(supplied by NCBI)'] and r['UniProt ID(supplied by UniProt)']:
                parts = r['UniProt ID(supplied by UniProt)'].split(';')
                for p in parts:
                    uniprot_to_entrez[p] = r['NCBI Gene ID(supplied by NCBI)']

    dropped_samples = [
        'LFQ intensity ARFGAP2_KO_3', 'LFQ intensity ARID1A_c010_KO_2',
        'LFQ intensity ARID1A_c010_KO_3', 'LFQ intensity ARID1B_KO_1',
        'LFQ intensity ARID1B_KO_2', 'LFQ intensity DDX5_KO_4'
     ] # Problem running these, technical replicates run in their place

    hap1_processed = process_proteomics_data(
        filename_input, dropped_samples=dropped_samples,
        entrez_map=uniprot_to_entrez, symbol_map=uniprot_to_symbol,
        drop_missing=False, fill=True, log_transform=True,
        percent_present=.2, adjust_means=True
    )
    lfq_columns = [x for x in hap1_processed.columns if 'LFQ' in x]
    hap1_short = hap1_processed[["GeneNames"] + lfq_columns].set_index("GeneNames")
    hap1_short = hap1_short.reset_index().rename(columns={"GeneNames": "gene_name"}).set_index("gene_name").rename(columns=lambda x: drop_lfq(x))

    HAP1_prot_renamed = hap1_short.rename(columns={'ARID1B_12': 'ARID1B_2'})

    all_clones = list(set(['_'.join(x.split('_')[0:-1]) for x in HAP1_prot_renamed.columns.to_list() if 'WT_' not in x]))
    all_clones = [x for x in all_clones if 'ASF1' not in x]

    A2_df = pd.DataFrame(all_clones).rename(columns={0: 'A2'})
    A2_df['t_stat'] = A2_df.A2.apply(lambda x: run_ttest_for_clone(x, df=HAP1_prot_renamed)[0])
    A2_df['p_val'] = A2_df.A2.apply(lambda x: run_ttest_for_clone(x, df=HAP1_prot_renamed)[1])
    A2_df = A2_df.dropna(subset='p_val')
    A2_df['p_values_adjusted'] = multipletests(A2_df['p_val'], method='fdr_bh')[1]
    A2_df['sig'] = (A2_df['p_val'] < 0.05) & (A2_df['t_stat'] < 0)

    data_A2 = {}
    for A2 in A2_df.A2.to_list():
        A2_noclone = A2.split('_')[0] # In case we have a clone ID appended to the gene name as we do when we have multiple viable clones
        genelist = [x for x in HAP1_prot_renamed.index.to_list() if ((f';{A2_noclone}' in x) or (f'{A2_noclone};' in x) or (x == A2_noclone))]
        gene_index = genelist[0]
        KO_samples = [x for x in HAP1_prot_renamed.columns.to_list() if A2 in x]
        KO_quants = HAP1_prot_renamed[(HAP1_prot_renamed.index == gene_index)][KO_samples].iloc[0].values.tolist()
        WT_quants = HAP1_prot_renamed[(HAP1_prot_renamed.index == gene_index)][[x for x in HAP1_prot_renamed.columns.to_list() if 'WT_' in x]].iloc[0].values.tolist()
        KO_quants = [x for x in KO_quants if not np.isnan(x)]
        WT_quants = [x for x in WT_quants if not np.isnan(x)]
        data_A2[A2] = KO_quants, WT_quants

    A2_df['logFC'] = A2_df['A2'].apply(lambda key: calculate_logFC(data_A2, key))
    A2_df_wclones = deepcopy(A2_df)
    A2_df = group_A2_clones(A2_df)
    A2_df['logFC'] = A2_df['logFC'].apply(lambda x: round(x, 3))
    valKOs = A2_df[A2_df.sig == True].A2.to_list()

    all_pairs = pd.read_csv('../data/general/all_pairs_wsexchr.csv', index_col=0)
    all_pairs = all_pairs[all_pairs.sex_chr == False] # Excluding sex chromosome paralogs- but none are testable in this analysis anyway 

    val_paralogs = {}
    for valKO in valKOs:
        val_paralogs[valKO] = all_pairs[all_pairs.A2 == (valKO.split('_')[0])].A1.to_list()

    pairs_to_test = [f'{x}_{k}' for k, v in val_paralogs.items() for x in v]

    A1A2_df = pd.DataFrame(pairs_to_test).rename(columns={0: 'gene_pair'})
    A1A2_df['t_stat'] = A1A2_df.gene_pair.apply(lambda x: run_ttest_for_paralog(x, df=HAP1_prot_renamed, A2_df_input=A2_df_wclones, A2_df_input_agg=A2_df)[0])
    A1A2_df['p_val'] = A1A2_df.gene_pair.apply(lambda x: run_ttest_for_paralog(x, df=HAP1_prot_renamed, A2_df_input=A2_df_wclones, A2_df_input_agg=A2_df)[1])
    A1A2_df = A1A2_df.dropna(subset='p_val')
    A1A2_df['p_values_adjusted'] = multipletests(A1A2_df['p_val'], method='fdr_bh')[1]
    A1A2_df['compensation'] = (A1A2_df['p_values_adjusted'] < 0.05) & (A1A2_df['p_val'] < 0.05) & (A1A2_df['t_stat'] > 0)
    A1A2_df['collateral_loss'] = (A1A2_df['p_values_adjusted'] < 0.05) & (A1A2_df['p_val'] < 0.05) & (A1A2_df['t_stat'] < 0)
    
    data = {}
    for pair in A1A2_df.gene_pair.to_list():
        splitlist = pair.split('_')
        A1 = splitlist[0]
        if len(splitlist) == 2:
            A2 = splitlist[1]
        elif len(splitlist) == 3:
            A2 = splitlist[1] + '_' + splitlist[2]
        genelist = [x for x in HAP1_prot_renamed.index.to_list() if ((f';{A1}' in x) or (f'{A1};' in x) or (x == A1))]
        if len(genelist)>0:
            gene_index = genelist[0]
        KO_samples = [x for x in HAP1_prot_renamed.columns.to_list() if A2 in x]
        KO_quants = HAP1_prot_renamed[(HAP1_prot_renamed.index == gene_index)][KO_samples].iloc[0].values.tolist()
        WT_quants = HAP1_prot_renamed[(HAP1_prot_renamed.index == gene_index)][[x for x in HAP1_prot_renamed.columns.to_list() if 'WT_' in x]].iloc[0].values.tolist()
        KO_quants = [x for x in KO_quants if np.isnan(x) == False]
        WT_quants = [x for x in WT_quants if np.isnan(x) == False]
        data[pair] = KO_quants, WT_quants
        A1A2_df['logFC'] = A1A2_df['gene_pair'].apply(lambda key: calculate_logFC(data, key))
        A1A2_df['logFC'] = A1A2_df['logFC'].round(3)

    A1A2_df['logFC'] = A1A2_df['gene_pair'].apply(lambda key: calculate_logFC(data, key))

    print(f'Final processed HAP1 proteomics dataset has {HAP1_prot_renamed.shape[0]} proteins and {HAP1_prot_renamed.shape[1]} unique samples incl WTs.')
    print(f'{len(A2_df)} A2s tested, {len(A2_df[A2_df.sig])} validate.')
    print(f'{len(A1A2_df)} paralog pairs tested, {len(A1A2_df[A1A2_df.compensation])} compensation and {len(A1A2_df[A1A2_df.collateral_loss])} collateral loss hits.')
    print(f'Compensation hits: (Gene pairs are in A1_A2 format where A2 is lost and A1 is the paralog)\n{",".join(A1A2_df[A1A2_df.compensation].gene_pair.to_list())}')
    print(f'Collateral loss hits:\n{",".join(A1A2_df[A1A2_df.collateral_loss].gene_pair.to_list())}')

    # Load overlap datasets
    dfs_for_cat_overlap, dfs_for_quant_overlap = load_and_process_data(overlap_only=True)

    # Generate overlaps for HAP1 results
    hap1_results = A1A2_df.copy()
    hap1_results['dataset'] = 'HAP1'
    hap1_results['FDR_threshold'] = 0.05
    hap1_results['gene_pair'] = hap1_results['gene_pair'].apply(lambda x: '_'.join(x.split('_')[0:2]))
    hap1_results['A1'] = hap1_results['gene_pair'].apply(lambda x: x.split('_')[0])
    hap1_results['A2'] = hap1_results['gene_pair'].apply(lambda x: x.split('_')[1])
    
    overlap_dfs = generate_overlaps(
        dfs_for_cat_overlap, 
        dfs_for_quant_overlap, 
        hap1_results, 
        runwith='prot', 
        workers=1, 
        skip_tests=True
    )

    # Unpack the overlap dataframes
    overlap_df_output, annotated_results = overlap_dfs

    # Save datasets
    A2_df.to_csv('../output/output_HAP1/HAP1_selftest_results.csv')
    A1A2_df.to_csv('../output/output_HAP1/HAP1_paralogtest_results.csv')

    HAP1_prot_renamed.index.name = 'gene_name'
    HAP1_prot_renamed.columns.name = ''
    HAP1_prot_renamed.to_csv('../output/output_HAP1/HAP1_prot_renamed.csv')

    # Save overlap dataframes
    overlap_df_output.to_csv('../output/output_HAP1/HAP1_overlaps_categorical.csv')
    annotated_results.to_csv('../output/output_HAP1/HAP1_tested_pairs_annots.csv')

if __name__ == "__main__":
    main()
