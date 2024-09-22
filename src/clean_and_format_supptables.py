# This script cleans output datasets and generates formatted supplementary tables

import numpy as np
import pandas as pd

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

    print('Formatting supplementary tables...')

    # Read and process HAP1 proteomics data
    supp_table1_HAP1_protdata = pd.read_csv('../output/output_HAP1/HAP1_prot_renamed.csv')
    supp_table1_HAP1_protdata = supp_table1_HAP1_protdata.rename(columns={'gene_name': 'gene_symbols'}).set_index('gene_symbols')
    supp_table1_HAP1_protdata.to_csv('../output/supp_tables/supplementary_table1_HAP1_proteomics.csv', index=True)

    # Process HAP1 self-test results
    supp_table2_HAP1_A2A2 = pd.read_csv('../output/output_HAP1/HAP1_selftest_results.csv', index_col=0)
    supp_table2_HAP1_A2A2 = supp_table2_HAP1_A2A2.rename(columns={'real_A2': 'gene_name', 'p_values_adjusted': 'FDR', 'sig': 'drop_in_KO'})
    supp_table2_HAP1_A2A2 = supp_table2_HAP1_A2A2.drop(columns=['A2', 'clone']).set_index('gene_name').reset_index()
    supp_table2_HAP1_A2A2.to_csv('../output/supp_tables/supplementary_table2_HAP1_selftest_results.csv', index=True)

    # Process HAP1 paralog test results
    supp_table3_HAP1_A1A2 = pd.read_csv('../output/output_HAP1/HAP1_paralogtest_results.csv', index_col=0)
    supp_table3_HAP1_A1A2 = supp_table3_HAP1_A1A2.rename(columns={'p_values_adjusted': 'FDR'})
    supp_table3_HAP1_A1A2['gene_pair'] = supp_table3_HAP1_A1A2['gene_pair'].apply(lambda x: '_'.join(x.split('_')[0:2]))
    supp_table3_HAP1_A1A2.to_csv('../output/supp_tables/supplementary_table3_HAP1_paralogtest_results.csv', index=True)

    # Process CPTAC self-test results
    supp_table4_CPTAC_A2A2 = pd.read_csv('../output/output_CPTAC/prot/self_tests_prot.csv', index_col=0)
    supp_table4_CPTAC_A2A2 = supp_table4_CPTAC_A2A2.rename(columns={'A2': 'A2_gene_symbol', 'backed': 'drop_when_lost', 'A2_lost_mean_quant_A1': 'A2_mean_when_lost', 'A2_other_mean_quant_A1': 'A2_mean_when_notlost'})
    supp_table4_CPTAC_A2A2.to_csv('../output/supp_tables/supplementary_table4_CPTAC_selftest_results.csv', index=True)

    # Process CPTAC proteomics paralog test results
    supp_table5_CPTAC_A1A2 = pd.read_csv('../output/output_CPTAC/prot/paralog_tests_prot.csv', index_col=0)
    supp_table5_CPTAC_A1A2 = supp_table5_CPTAC_A1A2.rename(columns={'p_adj': 'FDR'})
    supp_table5_CPTAC_A1A2['data_type'] = 'proteomics'
    supp_table5_CPTAC_A1A2.to_csv('../output/supp_tables/supplementary_table5_CPTAC_proteomics_paralogtest_results.csv', index=True)

    # Process transcriptomics and residual results
    trans_results = pd.read_csv('../output/output_CPTAC/trans/paralog_tests_trans.csv', index_col=0)
    trans_results = trans_results.rename(columns={'p_adj': 'FDR'})
    trans_results['data_type'] = 'transcriptomics'

    resid_results = pd.read_csv('../output/output_CPTAC/prot_residual/paralog_tests_prot_residual.csv', index_col=0)
    resid_results = resid_results.rename(columns={'p_adj': 'FDR'})
    resid_results['data_type'] = 'prot_residual'

    supp_table6_other2_datasets = pd.concat([trans_results, resid_results])
    supp_table6_other2_datasets['gene_pair_tested_in_dataset'] = supp_table6_other2_datasets['gene_pair'] + '_testedin_' + supp_table6_other2_datasets['data_type']
    supp_table6_other2_datasets = supp_table6_other2_datasets.set_index('gene_pair_tested_in_dataset').reset_index()
    supp_table6_other2_datasets.to_csv('../output/supp_tables/supplementary_table6_CPTAC_trans_and_resid_paralog_test_results.csv', index=True)

    # Process HAP1 biological annotations
    supp_table7_HAP1_bioannots = pd.read_csv('../output/output_HAP1/HAP1_overlaps_categorical.csv', index_col=0).set_index('gene_pair').iloc[:, 10:].reset_index()
    all_screened_pairs = pd.read_csv("../data/for_overlap/all_screened_paralog_pairs_25_04_22.csv")
    all_screened_pairs = all_screened_pairs[['sorted_gene_pair', 'n_SL_thompson', 'n_SL_dede', 'n_SL_parrish', 'n_SL_chymera', 'n_SL_ito']]

    renaming_dict = {
        'bronze_standard_SL': 'depmap_SL',
        'strict_comb_hit': 'SL_atleast2_CRISPRscreens',
        'in_PC_CORUM_essential': 'either_in_essential_CORUM_complex',
        'lenient_comb_hit': 'SL_atleast1_CRISPRscreen',
        'in_PC_CORUM': 'either_in_CORUM_complex',
        'in_PC_CORUM_both': 'both_in_same_CORUM_complex',
        'in_PC_EBI': 'either_in_EBIComplexPortal_complex',
        'in_PC_humap': 'either_in_humap_complex',
        'closest_pair': 'closest_paralogs_in_family',
        'famsize2': 'family_size2'
    }

    supp_table7_HAP1_bioannots = clean_up_overlap_annots(supp_table7_HAP1_bioannots, consensus_SL=all_screened_pairs, renaming_dict=renaming_dict)
    annotated_hap1 = pd.read_csv('../output/output_HAP1/HAP1_tested_pairs_annots.csv').drop(
        columns=['Unnamed: 0', 't_stat', 'p_val', 'p_values_adjusted', 'compensation', 'collateral_loss', 'logFC', 'dataset', 'FDR_threshold', 'sorted_gene_pair', 'A1', 'A2']).rename(columns={'category': 'HAP1_hit_type'})

    supp_table7_HAP1_bioannots = supp_table7_HAP1_bioannots.merge(annotated_hap1, on='gene_pair')
    supp_table7_HAP1_bioannots.to_csv('../output/supp_tables/supplementary_table7_allHAP1pairs_biological_info.csv')

    # Process CPTAC biological annotations
    supp_table8_CPTAC_bioannots = pd.read_csv('../output/output_CPTAC/prot/categorical_overlaps_prot.csv', index_col=0)
    annotated_cptac = pd.read_csv('../output/output_CPTAC/prot/all_quantoverlaps_prot.csv', index_col=0).drop(
        columns=['A1', 'A2', 'ols_p', 'ols_coef', 'n_A2_lost', 'lost_mean_quant', 'other_mean_quant', 'p_adj', 'compensation', 'collateral_loss', 'sorted_gene_pair']).rename(columns={'category': 'CPTAC_hit_type'})
    supp_table8_CPTAC_bioannots = clean_up_overlap_annots(supp_table8_CPTAC_bioannots, consensus_SL=all_screened_pairs, renaming_dict=renaming_dict)
    supp_table8_CPTAC_bioannots = supp_table8_CPTAC_bioannots.merge(annotated_cptac, on='gene_pair').set_index('gene_pair').reset_index()

    prot_tested = supp_table5_CPTAC_A1A2.gene_pair.to_list()
    protcomp = supp_table5_CPTAC_A1A2[supp_table5_CPTAC_A1A2.compensation].gene_pair.to_list()
    protcl = supp_table5_CPTAC_A1A2[supp_table5_CPTAC_A1A2.collateral_loss].gene_pair.to_list()

    trans_tested = trans_results.gene_pair.to_list()
    transcomp = trans_results[trans_results.compensation].gene_pair.to_list()
    transcl = trans_results[trans_results.collateral_loss].gene_pair.to_list()

    resid_tested = resid_results.gene_pair.to_list()
    residcomp = resid_results[resid_results.compensation].gene_pair.to_list()
    residcl = resid_results[resid_results.collateral_loss].gene_pair.to_list()

    supp_table8_CPTAC_bioannots['prot_compensation'] = supp_table8_CPTAC_bioannots['gene_pair'].apply(lambda x: (x in protcomp))
    supp_table8_CPTAC_bioannots['prot_collateral_loss'] = supp_table8_CPTAC_bioannots['gene_pair'].apply(lambda x: (x in protcl))
    supp_table8_CPTAC_bioannots['trans_compensation'] = supp_table8_CPTAC_bioannots['gene_pair'].apply(lambda x: (x in transcomp))
    supp_table8_CPTAC_bioannots['trans_collateral_loss'] = supp_table8_CPTAC_bioannots['gene_pair'].apply(lambda x: (x in transcl))
    supp_table8_CPTAC_bioannots['resid_compensation'] = supp_table8_CPTAC_bioannots['gene_pair'].apply(lambda x: (x in residcomp))
    supp_table8_CPTAC_bioannots['resid_collateral_loss'] = supp_table8_CPTAC_bioannots['gene_pair'].apply(lambda x: (x in residcl))

    supp_table8_CPTAC_bioannots['other_dir_prot_compensation'] = supp_table8_CPTAC_bioannots['gene_pair'].apply(lambda x: get_other_dir_genepair(x) in protcomp if get_other_dir_genepair(x) in prot_tested else np.nan)
    supp_table8_CPTAC_bioannots['other_dir_prot_collateral_loss'] = supp_table8_CPTAC_bioannots['gene_pair'].apply(lambda x: get_other_dir_genepair(x) in protcl if get_other_dir_genepair(x) in prot_tested else np.nan)
    supp_table8_CPTAC_bioannots['other_dir_trans_compensation'] = supp_table8_CPTAC_bioannots['gene_pair'].apply(lambda x: get_other_dir_genepair(x) in transcomp if get_other_dir_genepair(x) in trans_tested else np.nan)
    supp_table8_CPTAC_bioannots['other_dir_trans_collateral_loss'] = supp_table8_CPTAC_bioannots['gene_pair'].apply(lambda x: get_other_dir_genepair(x) in transcl if get_other_dir_genepair(x) in trans_tested else np.nan)
    supp_table8_CPTAC_bioannots['other_dir_resid_compensation'] = supp_table8_CPTAC_bioannots['gene_pair'].apply(lambda x: get_other_dir_genepair(x) in residcomp if get_other_dir_genepair(x) in resid_tested else np.nan)
    supp_table8_CPTAC_bioannots['other_dir_resid_collateral_loss'] = supp_table8_CPTAC_bioannots['gene_pair'].apply(lambda x: get_other_dir_genepair(x) in residcl if get_other_dir_genepair(x) in resid_tested else np.nan)

    supp_table8_CPTAC_bioannots.to_csv('../output/supp_tables/supplementary_table8_allCPTACpairs_biological_info.csv')

    # Process CPTAC Fisher's Exact Test results
    supp_table9_CPTAC_FET = pd.read_csv('../output/output_CPTAC/prot/categorical_FETs_uniquegenepairs_prot.csv', index_col=0)
    supp_table9_CPTAC_FET['interaction_dataset'] = supp_table9_CPTAC_FET['interaction_dataset'].apply(lambda x: renaming_dict[x] if x in renaming_dict.keys() else x)
    supp_table9_CPTAC_FET['dataset'] = 'proteomics'

    trans_FET = pd.read_csv('../output/output_CPTAC/trans/categorical_FETs_uniquegenepairs_trans.csv', index_col=0)
    trans_FET['interaction_dataset'] = trans_FET['interaction_dataset'].apply(lambda x: renaming_dict[x] if x in renaming_dict.keys() else x)
    trans_FET['dataset'] = 'transcriptomics'

    resid_FET = pd.read_csv('../output/output_CPTAC/prot_residual/categorical_FETs_uniquegenepairs_prot_residual.csv', index_col=0)
    resid_FET['interaction_dataset'] = resid_FET['interaction_dataset'].apply(lambda x: renaming_dict[x] if x in renaming_dict.keys() else x)
    resid_FET['dataset'] = 'prot_residual'

    supp_table9_CPTAC_FET = pd.concat([supp_table9_CPTAC_FET, trans_FET, resid_FET])
    supp_table9_CPTAC_FET.to_csv('../output/supp_tables/supplementary_table9_allCPTAC_categorical_overlaptests.csv', index=True)

    # Process CPTAC t-test results
    supp_table10_CPTAC_tt = pd.read_csv('../output/output_CPTAC/prot/quantitative_ttests_foroverlaps_prot.csv', index_col=0)
    supp_table10_CPTAC_tt['dataset'] = 'proteomics'
    supp_table10_CPTAC_tt = supp_table10_CPTAC_tt.drop(columns='colname').set_index('Variable').reset_index().rename(columns={'Variable': 'variable'})
    supp_table10_CPTAC_tt.to_csv('../output/supp_tables/supplementary_table10_allCPTAC_quantitative_overlaptests.csv', index=True)

    print('Done!')

if __name__ == '__main__':
    main()