# This script cleans output datasets and generates formatted appendix tables in a multi-sheet Excel file

import numpy as np
import pandas as pd
from openpyxl import Workbook
from openpyxl.utils.dataframe import dataframe_to_rows
import os

def main():

    def get_other_dir_genepair(gene_pair):
        A1, A2 = gene_pair.split('_')
        return f'{A2}_{A1}'

    def clean_up_overlap_annots(df, consensus_SL, renaming_dict):
        df = df[list(renaming_dict.keys()) + ['gene_pair', 'sorted_gene_pair']]
        df = df.rename(columns=renaming_dict)
        consensus_SL = consensus_SL[['sorted_gene_pair', 'n_SL_thompson', 'n_SL_dede', 'n_SL_parrish', 'n_SL_chymera', 'n_SL_ito']]
        df = df.merge(consensus_SL, on='sorted_gene_pair', how='left').set_index('gene_pair').reset_index().set_index('sorted_gene_pair').reset_index()
        return df

    print('Formatting appendix tables...')

    # Create EV_datasets directory if it doesn't exist
    os.makedirs('../output/EV_datasets', exist_ok=True)
    
    # Define table descriptions for the introduction sheet
    table_descriptions = [
        "**Dataset EV1: Processed HAP1 proteomic data** HAP1 proteomic data from 34 paralog knockouts, processed as outlined in the Methods.",
        "**Dataset EV2: All self-abundance HAP1 tests** T-statistic, log fold change, p-values and FDRs are for a comparison of protein abundance when the gene has been knocked out, versus its abundance in the wild-type. \"drop_in_KO\" contains information on whether the gene passed the test.",
        "**Dataset EV3: All paralog HAP1 tests** T-statistic, log fold change, p-values and FDRs are for a comparison of protein abundance when its paralog has been knocked out, versus its abundance in the wild-type. \"compensation\" and \"collateral_loss\" contain information on whether each paralog is a compensation or collateral loss hit (see Methods).",
        "**Dataset EV4: All self-abundance CPTAC tests**: Ordinary least squares regression models were fit for each testable protein to explain its abundance using its own hemizygous loss status across CPTAC samples, with lineage/study as a covariate. OLS coefficients, p-values (i.e. Two-tailed p values for the t-statistic for the A2 loss variable), FDRs (Benjamini Hochberg multiple testing correction applied to p-values), mean values when lost vs. when not lost, number of samples a gene has been lost in, and whether or not loss is significantly (p-value < 0.05) associated with drop in protein abundance.",
        "**Dataset EV5: CPTAC proteomic results** All paralog CPTAC tests. Ordinary least squares regression models were fit for each testable protein to explain its abundance using its paralog's hemizygous loss status across CPTAC samples, with lineage/study as a covariate. OLS coefficients, p-values (i.e. t-statistics for the loss variable), FDRs (Benjamini Hochberg multiple testing correction applied to pvalues), mean values when lost vs. when not lost, and OLS model r-squared values, number of samples a gene has been lost in, and whether or not paralog loss is significantly (FDR < 5% and uncorrected p-value < 0.05) associated with change in protein abundance (compensation or collateral_loss).",
        "**Dataset EV6: CPTAC transcriptomic and protein residual results** Results of the above analysis run using transcriptomic data rather than proteomic data, as well as results with protein residuals, i.e. a version of the proteomic dataset where lineage and self-transcript effects have been regressed out (by fitting separate ordinary least squares models for each protein, see Methods).",
        "**Dataset EV7: Biological information for all HAP1 pairs** Information about protein complex membership, closest pair status, sequence identity, family size, Jaccard index, degree centrality, conservation score, and synthetic lethality for all paralog pairs tested using HAP1 data.",
        "**Dataset EV8: Biological information for all CPTAC pairs** Information about protein complex membership, closest pair status, sequence identity, family size, Jaccard index, degree centrality, conservation score, and synthetic lethality for all paralog pairs tested using CPTAC data.",
        "**Dataset EV9: All Fishers Exact Test results (categorical overlap tests) for CPTAC pairs** Results of all Fishers Exact Tests run to identify overlap between CPTAC compensation and collateral loss status and various categorical biological characteristics including synthetic lethality, protein complex membership, closest pair status, and family size 2 (i.e no other paralogs in the family) as described in the methods.",
        "**Dataset EV10: All t-test results (quantitative overlap tests) for CPTAC pairs** Results of t-tests run to identify enrichment in CPTAC hits for quantitative biological characteristics such as sequence identity, family size, Jaccard index, degree centrality of lost gene, essentiality of the neighbours of the lost gene, and conservation scores."
    ]

    # Read and process HAP1 proteomics data
    appendix_table1_HAP1_protdata = pd.read_csv('../output/output_HAP1/HAP1_prot_renamed.csv')
    appendix_table1_HAP1_protdata = appendix_table1_HAP1_protdata.rename(columns={'gene_name': 'gene_symbols'}).set_index('gene_symbols')
    appendix_table1_HAP1_protdata.to_csv('../output/appendix_tables/appendix_table1_HAP1_proteomics.csv', index=True)

    # Process HAP1 self-test results
    appendix_table2_HAP1_A2A2 = pd.read_csv('../output/output_HAP1/HAP1_selftest_results.csv', index_col=0)
    appendix_table2_HAP1_A2A2 = appendix_table2_HAP1_A2A2.rename(columns={'real_A2': 'gene_name', 'p_values_adjusted': 'FDR', 'sig': 'drop_in_KO'})
    appendix_table2_HAP1_A2A2 = appendix_table2_HAP1_A2A2.drop(columns=['A2', 'clone']).set_index('gene_name').reset_index()
    # Ensure boolean columns stay as strings when saving
    if 'drop_in_KO' in appendix_table2_HAP1_A2A2.columns:
        appendix_table2_HAP1_A2A2['drop_in_KO'] = appendix_table2_HAP1_A2A2['drop_in_KO'].map({True: 'True', False: 'False'})
    appendix_table2_HAP1_A2A2.to_csv('../output/appendix_tables/appendix_table2_HAP1_selftest_results.csv', index=True)

    # Process HAP1 paralog test results
    appendix_table3_HAP1_A1A2 = pd.read_csv('../output/output_HAP1/HAP1_paralogtest_results.csv', index_col=0)
    appendix_table3_HAP1_A1A2 = appendix_table3_HAP1_A1A2.rename(columns={'p_values_adjusted': 'FDR'})
    appendix_table3_HAP1_A1A2['gene_pair'] = appendix_table3_HAP1_A1A2['gene_pair'].apply(lambda x: '_'.join(x.split('_')[0:2]))
    # Ensure boolean columns stay as strings when saving
    for col in ['compensation', 'collateral_loss']:
        if col in appendix_table3_HAP1_A1A2.columns:
            appendix_table3_HAP1_A1A2[col] = appendix_table3_HAP1_A1A2[col].map({True: 'True', False: 'False'})
    appendix_table3_HAP1_A1A2.to_csv('../output/appendix_tables/appendix_table3_HAP1_paralogtest_results.csv', index=True)

    # Process CPTAC self-test results
    appendix_table4_CPTAC_A2A2 = pd.read_csv('../output/output_CPTAC/prot/self_tests_prot.csv', index_col=0)
    appendix_table4_CPTAC_A2A2 = appendix_table4_CPTAC_A2A2.rename(columns={'A2': 'A2_gene_symbol', 'backed': 'drop_when_lost', 'A2_lost_mean_quant_A1': 'A2_mean_when_lost', 'A2_other_mean_quant_A1': 'A2_mean_when_notlost'})
    # Ensure boolean columns stay as strings when saving
    if 'drop_when_lost' in appendix_table4_CPTAC_A2A2.columns:
        appendix_table4_CPTAC_A2A2['drop_when_lost'] = appendix_table4_CPTAC_A2A2['drop_when_lost'].map({True: 'True', False: 'False'})
    appendix_table4_CPTAC_A2A2.to_csv('../output/appendix_tables/appendix_table4_CPTAC_selftest_results.csv', index=True)

    # Process CPTAC proteomics paralog test results
    appendix_table5_CPTAC_A1A2 = pd.read_csv('../output/output_CPTAC/prot/paralog_tests_prot.csv', index_col=0)
    appendix_table5_CPTAC_A1A2 = appendix_table5_CPTAC_A1A2.rename(columns={'p_adj': 'FDR'})
    appendix_table5_CPTAC_A1A2['data_type'] = 'proteomics'
    # Ensure boolean columns stay as strings when saving
    for col in ['compensation', 'collateral_loss']:
        if col in appendix_table5_CPTAC_A1A2.columns:
            appendix_table5_CPTAC_A1A2[col] = appendix_table5_CPTAC_A1A2[col].map({True: 'True', False: 'False'})
    appendix_table5_CPTAC_A1A2.to_csv('../output/appendix_tables/appendix_table5_CPTAC_proteomics_paralogtest_results.csv', index=True)

    # Process transcriptomics and residual results
    trans_results = pd.read_csv('../output/output_CPTAC/trans/paralog_tests_trans.csv', index_col=0)
    trans_results = trans_results.rename(columns={'p_adj': 'FDR'})
    trans_results['data_type'] = 'transcriptomics'
    # Ensure boolean columns stay as strings
    for col in ['compensation', 'collateral_loss']:
        if col in trans_results.columns:
            trans_results[col] = trans_results[col].map({True: 'True', False: 'False'})

    resid_results = pd.read_csv('../output/output_CPTAC/prot_residual/paralog_tests_prot_residual.csv', index_col=0)
    resid_results = resid_results.rename(columns={'p_adj': 'FDR'})
    resid_results['data_type'] = 'prot_residual'
    # Ensure boolean columns stay as strings
    for col in ['compensation', 'collateral_loss']:
        if col in resid_results.columns:
            resid_results[col] = resid_results[col].map({True: 'True', False: 'False'})

    appendix_table6_other2_datasets = pd.concat([trans_results, resid_results])
    appendix_table6_other2_datasets['gene_pair_tested_in_dataset'] = appendix_table6_other2_datasets['gene_pair'] + '_testedin_' + appendix_table6_other2_datasets['data_type']
    appendix_table6_other2_datasets = appendix_table6_other2_datasets.set_index('gene_pair_tested_in_dataset').reset_index()
    appendix_table6_other2_datasets.to_csv('../output/appendix_tables/appendix_table6_CPTAC_trans_and_resid_paralog_test_results.csv', index=True)

    # Process HAP1 biological annotations
    appendix_table7_HAP1_bioannots = pd.read_csv('../output/output_HAP1/HAP1_overlaps_categorical.csv', index_col=0).set_index('gene_pair').iloc[:, 10:].reset_index()
    all_screened_pairs = pd.read_csv("../data/for_overlap/all_screened_paralog_pairs_25_04_22.csv")
    all_screened_pairs = all_screened_pairs[['sorted_gene_pair', 'n_SL_thompson', 'n_SL_dede', 'n_SL_parrish', 'n_SL_chymera', 'n_SL_ito']]

    renaming_dict = {
        'bronze_standard_SL': 'depmap_SL',
        'strict_comb_hit': 'SL_atleast2_CRISPRscreens',
        'lenient_comb_hit': 'SL_atleast1_CRISPRscreen',
        'in_PC_CORUM_essential': 'either_in_essential_CORUM_complex',
        'in_PC_CORUM': 'either_in_CORUM_complex',
        'in_PC_CORUM_both': 'both_in_same_CORUM_complex',
        'in_PC_EBI': 'either_in_EBIComplexPortal_complex',
        'in_PC_humap': 'either_in_humap_complex',
        'closest_pair': 'closest_paralogs_in_family',
        'famsize2': 'family_size2'
    }

    appendix_table7_HAP1_bioannots = clean_up_overlap_annots(appendix_table7_HAP1_bioannots, consensus_SL=all_screened_pairs, renaming_dict=renaming_dict)
    annotated_hap1 = pd.read_csv('../output/output_HAP1/HAP1_tested_pairs_annots.csv').drop(
        columns=['Unnamed: 0', 't_stat', 'p_val', 'p_values_adjusted', 'compensation', 'collateral_loss', 'logFC', 'dataset', 'FDR_threshold', 'sorted_gene_pair', 'A1', 'A2']).rename(columns={'category': 'HAP1_hit_type'})

    appendix_table7_HAP1_bioannots = appendix_table7_HAP1_bioannots.merge(annotated_hap1, on='gene_pair')
    
    # Convert boolean columns to string True/False
    for col in appendix_table7_HAP1_bioannots.columns:
        if appendix_table7_HAP1_bioannots[col].dtype == bool:
            appendix_table7_HAP1_bioannots[col] = appendix_table7_HAP1_bioannots[col].map({True: 'True', False: 'False'})
    
    appendix_table7_HAP1_bioannots.to_csv('../output/appendix_tables/appendix_table7_allHAP1pairs_biological_info.csv')

    # Process CPTAC biological annotations
    appendix_table8_CPTAC_bioannots = pd.read_csv('../output/output_CPTAC/prot/categorical_overlaps_prot.csv', index_col=0)
    annotated_cptac = pd.read_csv('../output/output_CPTAC/prot/all_quantoverlaps_prot.csv', index_col=0).drop(
        columns=['A1', 'A2', 'ols_p', 'ols_coef', 'n_A2_lost', 'lost_mean_quant', 'other_mean_quant', 'p_adj', 'compensation', 'collateral_loss', 'sorted_gene_pair']).rename(columns={'category': 'CPTAC_hit_type'})
    appendix_table8_CPTAC_bioannots = clean_up_overlap_annots(appendix_table8_CPTAC_bioannots, consensus_SL=all_screened_pairs, renaming_dict=renaming_dict)
    appendix_table8_CPTAC_bioannots = appendix_table8_CPTAC_bioannots.merge(annotated_cptac, on='gene_pair').set_index('gene_pair').reset_index()

    # Need to use string versions of boolean values to avoid conversion to 0/1
    prot_tested = appendix_table5_CPTAC_A1A2.gene_pair.to_list()
    protcomp = appendix_table5_CPTAC_A1A2[appendix_table5_CPTAC_A1A2.compensation == 'True'].gene_pair.to_list()
    protcl = appendix_table5_CPTAC_A1A2[appendix_table5_CPTAC_A1A2.collateral_loss == 'True'].gene_pair.to_list()

    trans_tested = trans_results.gene_pair.to_list()
    transcomp = trans_results[trans_results.compensation == 'True'].gene_pair.to_list()
    transcl = trans_results[trans_results.collateral_loss == 'True'].gene_pair.to_list()

    resid_tested = resid_results.gene_pair.to_list()
    residcomp = resid_results[resid_results.compensation == 'True'].gene_pair.to_list()
    residcl = resid_results[resid_results.collateral_loss == 'True'].gene_pair.to_list()

    # Create boolean columns but immediately convert to string True/False
    appendix_table8_CPTAC_bioannots['prot_compensation'] = appendix_table8_CPTAC_bioannots['gene_pair'].apply(lambda x: str(x in protcomp))
    appendix_table8_CPTAC_bioannots['prot_collateral_loss'] = appendix_table8_CPTAC_bioannots['gene_pair'].apply(lambda x: str(x in protcl))
    appendix_table8_CPTAC_bioannots['trans_compensation'] = appendix_table8_CPTAC_bioannots['gene_pair'].apply(lambda x: str(x in transcomp))
    appendix_table8_CPTAC_bioannots['trans_collateral_loss'] = appendix_table8_CPTAC_bioannots['gene_pair'].apply(lambda x: str(x in transcl))
    appendix_table8_CPTAC_bioannots['resid_compensation'] = appendix_table8_CPTAC_bioannots['gene_pair'].apply(lambda x: str(x in residcomp))
    appendix_table8_CPTAC_bioannots['resid_collateral_loss'] = appendix_table8_CPTAC_bioannots['gene_pair'].apply(lambda x: str(x in residcl))

    # Handle NaN values properly for 'other_dir' columns
    appendix_table8_CPTAC_bioannots['other_dir_prot_compensation'] = appendix_table8_CPTAC_bioannots['gene_pair'].apply(
        lambda x: str(get_other_dir_genepair(x) in protcomp) if get_other_dir_genepair(x) in prot_tested else np.nan)
    appendix_table8_CPTAC_bioannots['other_dir_prot_collateral_loss'] = appendix_table8_CPTAC_bioannots['gene_pair'].apply(
        lambda x: str(get_other_dir_genepair(x) in protcl) if get_other_dir_genepair(x) in prot_tested else np.nan)
    appendix_table8_CPTAC_bioannots['other_dir_trans_compensation'] = appendix_table8_CPTAC_bioannots['gene_pair'].apply(
        lambda x: str(get_other_dir_genepair(x) in transcomp) if get_other_dir_genepair(x) in trans_tested else np.nan)
    appendix_table8_CPTAC_bioannots['other_dir_trans_collateral_loss'] = appendix_table8_CPTAC_bioannots['gene_pair'].apply(
        lambda x: str(get_other_dir_genepair(x) in transcl) if get_other_dir_genepair(x) in trans_tested else np.nan)
    appendix_table8_CPTAC_bioannots['other_dir_resid_compensation'] = appendix_table8_CPTAC_bioannots['gene_pair'].apply(
        lambda x: str(get_other_dir_genepair(x) in residcomp) if get_other_dir_genepair(x) in resid_tested else np.nan)
    appendix_table8_CPTAC_bioannots['other_dir_resid_collateral_loss'] = appendix_table8_CPTAC_bioannots['gene_pair'].apply(
        lambda x: str(get_other_dir_genepair(x) in residcl) if get_other_dir_genepair(x) in resid_tested else np.nan)

    # Convert any remaining boolean columns to string True/False
    for col in appendix_table8_CPTAC_bioannots.columns:
        if appendix_table8_CPTAC_bioannots[col].dtype == bool:
            appendix_table8_CPTAC_bioannots[col] = appendix_table8_CPTAC_bioannots[col].map({True: 'True', False: 'False'})
            
    appendix_table8_CPTAC_bioannots.to_csv('../output/appendix_tables/appendix_table8_allCPTACpairs_biological_info.csv')

    # Process CPTAC Fisher's Exact Test results
    appendix_table9_CPTAC_FET = pd.read_csv('../output/output_CPTAC/prot/categorical_FETs_uniquegenepairs_prot.csv', index_col=0)
    appendix_table9_CPTAC_FET['interaction_dataset'] = appendix_table9_CPTAC_FET['interaction_dataset'].apply(lambda x: renaming_dict[x] if x in renaming_dict.keys() else x)
    appendix_table9_CPTAC_FET['dataset'] = 'proteomics'

    trans_FET = pd.read_csv('../output/output_CPTAC/trans/categorical_FETs_uniquegenepairs_trans.csv', index_col=0)
    trans_FET['interaction_dataset'] = trans_FET['interaction_dataset'].apply(lambda x: renaming_dict[x] if x in renaming_dict.keys() else x)
    trans_FET['dataset'] = 'transcriptomics'

    resid_FET = pd.read_csv('../output/output_CPTAC/prot_residual/categorical_FETs_uniquegenepairs_prot_residual.csv', index_col=0)
    resid_FET['interaction_dataset'] = resid_FET['interaction_dataset'].apply(lambda x: renaming_dict[x] if x in renaming_dict.keys() else x)
    resid_FET['dataset'] = 'prot_residual'

    appendix_table9_CPTAC_FET = pd.concat([appendix_table9_CPTAC_FET, trans_FET, resid_FET])
    appendix_table9_CPTAC_FET.to_csv('../output/appendix_tables/appendix_table9_allCPTAC_categorical_overlaptests.csv', index=True)

    # Process CPTAC t-test results
    appendix_table10_CPTAC_tt = pd.read_csv('../output/output_CPTAC/prot/quantitative_ttests_foroverlaps_prot.csv', index_col=0)
    appendix_table10_CPTAC_tt['dataset'] = 'proteomics'
    appendix_table10_CPTAC_tt = appendix_table10_CPTAC_tt.drop(columns='colname').set_index('Variable').reset_index().rename(columns={'Variable': 'variable'})
    appendix_table10_CPTAC_tt.to_csv('../output/appendix_tables/appendix_table10_allCPTAC_quantitative_overlaptests.csv', index=True)
    
    # Create individual Excel files for each dataset
    print('Creating individual Excel files for each dataset...')
    
    # Load all tables
    tables = [
        appendix_table1_HAP1_protdata.reset_index(),
        appendix_table2_HAP1_A2A2,
        appendix_table3_HAP1_A1A2,
        appendix_table4_CPTAC_A2A2,
        appendix_table5_CPTAC_A1A2,
        appendix_table6_other2_datasets,
        appendix_table7_HAP1_bioannots,
        appendix_table8_CPTAC_bioannots,
        appendix_table9_CPTAC_FET,
        appendix_table10_CPTAC_tt
    ]
    
    # Create and save individual files
    file_names = [
        "Dataset_EV1",
        "Dataset_EV2",
        "Dataset_EV3",
        "Dataset_EV4",
        "Dataset_EV5",
        "Dataset_EV6", 
        "Dataset_EV7",
        "Dataset_EV8",
        "Dataset_EV9",
        "Dataset_EV10"
    ]
    
    # Add tables to individual files
    for i, (table, file_name, description) in enumerate(zip(tables, file_names, table_descriptions)):
        # Create new workbook
        wb = Workbook()
        
        # Create legend sheet
        legend_sheet = wb.active
        legend_sheet.title = "Legend"
        legend_sheet.cell(row=1, column=1, value=description)
        
        # Create data sheet
        data_sheet = wb.create_sheet(title="Data")
        
        # Convert boolean columns to strings
        for col in table.columns:
            if table[col].dtype == bool:
                table[col] = table[col].map({True: 'True', False: 'False'})
                
        # Add data to data sheet
        for r_idx, row in enumerate(dataframe_to_rows(table, index=False, header=True), 1):
            for c_idx, value in enumerate(row, 1):
                if isinstance(value, bool):
                    value = str(value)
                data_sheet.cell(row=r_idx, column=c_idx, value=value)
        
        # Save the workbook
        wb.save(f'../output/EV_datasets/{file_name}.xlsx')
    
    print('Done! All datasets saved to ../output/EV_datasets/')

if __name__ == '__main__':
    main()
