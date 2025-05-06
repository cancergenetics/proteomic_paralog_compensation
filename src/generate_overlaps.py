import numpy as np 
import pandas as pd
import os
import scipy.stats as stats
from statsmodels.stats.multitest import fdrcorrection
from copy import deepcopy
import networkx as nx
from pandarallel import pandarallel
import requests
import time
from tqdm import tqdm
tqdm.pandas()        

def check_if_both_proteins_in_same_complex(input_gene_pair, corum_with_symbols):
    '''
    Checks if both proteins in the gene pair belong to the same CORUM complex
    '''
    A1, A2 = input_gene_pair.split("_")
    all_complexes_A1_list = corum_with_symbols[corum_with_symbols.gene_symbol == A1].complex.to_list()
    all_complexes_A1_df = corum_with_symbols[corum_with_symbols.complex.isin(all_complexes_A1_list)]
    all_complex_partners_A1 = all_complexes_A1_df.gene_symbol.to_list()
    return A2 in all_complex_partners_A1


def identify_essential_genes(df, threshold=-0.6, min_percentage=80):
    '''
    Filters DepMap processed CERES scores dataset to identify broadly essential genes 
    at a threshold of -0.6 
    '''
    total_cell_lines = len(df)
    min_count = int((min_percentage / 100) * total_cell_lines)
    essential_genes = []
    for gene in df.columns:
        if not pd.api.types.is_numeric_dtype(df[gene]):
            continue
        if (df[gene] < threshold).sum() >= min_count:
            essential_genes.append(gene)        
    return essential_genes


def run_FET(interaction_dataset, input_hit_type, gene_pair, overlap_df_output):
    '''
    Uses a generated overlap DF to run a Fishers Exact Test to check enrichment for a biological feature in 
    compensation and collateral loss hits (as compared to non hits)
    '''
    both_hit_types = {'compensation':'collateral_loss', 'collateral_loss':'compensation'}
    other_hit_type = both_hit_types[input_hit_type]
    overlap_df_local = deepcopy(overlap_df_output)
    overlap_df_local = overlap_df_local[
        (overlap_df_local[other_hit_type] == False)
    ][[gene_pair, interaction_dataset, input_hit_type]]
    overlap_df_local = overlap_df_local.dropna(subset=[input_hit_type, interaction_dataset]).drop_duplicates(subset=[gene_pair]).set_index(gene_pair)
    num_pairs_overlap = overlap_df_local.shape[0]
    crosstab2 = pd.crosstab(overlap_df_local[input_hit_type], overlap_df_local[interaction_dataset])
    if crosstab2.shape == (2,2):
        FET_OR_local = stats.fisher_exact(crosstab2)[0]
        FET_p_local = stats.fisher_exact(crosstab2)[1]
        
    else:
        FET_OR_local, FET_p_local = np.nan, np.nan
    return pd.DataFrame({
        'interaction_dataset': [interaction_dataset],
        'hit_type': [input_hit_type],
        'common_pairs': [num_pairs_overlap],
        'FET_OR': [FET_OR_local],
        'FET_p': [FET_p_local],
    })

def process_FET(inputdf):
    '''
    Corrects p-values for multiple testing to obtain FDR
    '''
    inputdf = inputdf.dropna(subset='FET_p')
    inputdf = inputdf.replace([np.inf, -np.inf], np.nan).dropna(subset=['FET_OR'])
    inputdf['FET_p_adj'] = fdrcorrection(inputdf.FET_p.values)[1]
    return inputdf

def build_ppi_network(df):
    '''
    Given a PPI dataframe, builds a network representation using NetworkX
    '''
    G = nx.Graph()
    for idx, row in tqdm(df.iterrows()):
        G.add_edge(row['gene_A'], row['gene_B'])
    return G

def jaccard_index(G, gene_pair):
    '''
    Given a PPI network represented as a graph, calculates the Jaccard index for a pair of nodes i.e. proteins
    '''
    protein1, protein2 = gene_pair.split('_')
    if not G.has_node(protein1) or not G.has_node(protein2):
        return np.nan  
    neighbors1 = set(G.neighbors(protein1))
    neighbors2 = set(G.neighbors(protein2))
    intersection = neighbors1.intersection(neighbors2)
    union = neighbors1.union(neighbors2)
    if not union:
        return 0
    return len(intersection) / len(union)

def get_ppi_ess(gene, network, scores):
    '''
    Given a node i.e. protein in the PPI network, returns median essentiality (CERES score) of its immediate neighbours
    '''
    if gene not in network.nodes:
        return np.nan  
    neighbours = list(set(network.neighbors(gene)))
    all_neighbour_scores = []
    for genes in neighbours:
        if genes in scores.columns:
            all_neighbour_scores.extend(scores[genes].tolist()) 
    return np.nanmedian(all_neighbour_scores) if all_neighbour_scores else np.nan

def multiple_ttests_with_correction(df, colnames, drop):
    '''
    Runs multiple t-tests to identify enrichment for quantitative biological variables,
    applies multiple testing correction, includes degrees of freedom and sample size for each test,
    and calculates Cohen's d effect size
    '''
    results = []
    
    for colname in colnames:
        if drop:
            df2 = df.drop_duplicates(subset='sorted_gene_pair')
        else:
            df2 = df
        compvals = df2[df2.compensation].dropna(subset=colname)[colname].to_list()
        clvals = df2[df2.collateral_loss].dropna(subset=colname)[colname].to_list()
        nonvals = df2[df2.category == 'non_hit'].dropna(subset=colname)[colname].to_list()
        
        compvals = [x for x in compvals if not np.isinf(x)]
        clvals = [x for x in clvals if not np.isinf(x)]
        nonvals = [x for x in nonvals if not np.isinf(x)]

        comp_non = stats.ttest_ind(compvals, nonvals)
        cl_non = stats.ttest_ind(clvals, nonvals)
        
        comp_non_df = len(compvals) + len(nonvals) - 2
        cl_non_df = len(clvals) + len(nonvals) - 2
        comp_non_sample_size = len(compvals) + len(nonvals)
        cl_non_sample_size = len(clvals) + len(nonvals)
        
        # Calculate Cohen's d for compensation vs non-hit
        if len(compvals) > 0 and len(nonvals) > 0:
            # Calculate means
            comp_mean = np.mean(compvals)
            non_mean = np.mean(nonvals)
            
            # Calculate pooled standard deviation
            comp_var = np.var(compvals, ddof=1)  # ddof=1 for sample variance
            non_var = np.var(nonvals, ddof=1)
            
            # Pooled SD formula: sqrt(((n1-1)*s1^2 + (n2-1)*s2^2) / (n1+n2-2))
            pooled_sd_comp_non = np.sqrt(
                ((len(compvals) - 1) * comp_var + (len(nonvals) - 1) * non_var) / 
                (len(compvals) + len(nonvals) - 2)
            )
            
            # Cohen's d = (mean1 - mean2) / pooled_sd
            cohens_d_comp_non = (comp_mean - non_mean) / pooled_sd_comp_non if pooled_sd_comp_non > 0 else 0
        else:
            cohens_d_comp_non = np.nan
        
        # Calculate Cohen's d for collateral loss vs non-hit
        if len(clvals) > 0 and len(nonvals) > 0:
            # Calculate means
            cl_mean = np.mean(clvals)
            non_mean = np.mean(nonvals)
            
            # Calculate pooled standard deviation
            cl_var = np.var(clvals, ddof=1)
            non_var = np.var(nonvals, ddof=1)
            
            # Pooled SD formula
            pooled_sd_cl_non = np.sqrt(
                ((len(clvals) - 1) * cl_var + (len(nonvals) - 1) * non_var) / 
                (len(clvals) + len(nonvals) - 2)
            )
            
            # Cohen's d
            cohens_d_cl_non = (cl_mean - non_mean) / pooled_sd_cl_non if pooled_sd_cl_non > 0 else 0
        else:
            cohens_d_cl_non = np.nan
        
        results.append({
            'colname': colname,
            'n_comp_values': len(compvals),
            'n_non_values': len(nonvals),
            'n_cl_values': len(clvals),
            'comp_vs_non_tstat': comp_non.statistic,
            'comp_vs_non_pvalue': comp_non.pvalue,
            'comp_vs_non_df': comp_non_df,
            'comp_vs_non_sample_size': comp_non_sample_size,
            'comp_vs_non_cohens_d': cohens_d_comp_non,  # Added Cohen's d
            'cl_vs_non_tstat': cl_non.statistic,
            'cl_vs_non_pvalue': cl_non.pvalue,
            'cl_vs_non_df': cl_non_df,
            'cl_vs_non_sample_size': cl_non_sample_size,
            'cl_vs_non_cohens_d': cohens_d_cl_non  # Added Cohen's d
        })
    
    results_df = pd.DataFrame(results)
    
    _, comp_corrected_pvalues = fdrcorrection(results_df['comp_vs_non_pvalue'])
    _, cl_corrected_pvalues = fdrcorrection(results_df['cl_vs_non_pvalue'])
    
    results_df['comp_vs_non_FDR'] = comp_corrected_pvalues
    results_df['cl_vs_non_FDR'] = cl_corrected_pvalues
    
    return results_df

def get_hit_type(comp, cl):
    if comp:
        return 'compensation'
    elif cl:
        return 'collateral_loss'
    else:
        return 'non_hit'

def generate_overlaps(dfs_for_cat_overlap, dfs_for_quant_overlap, results, runwith, workers, skip_tests=False):
    '''
    Given a results dataframe with a set of compensation and collateral loss pairs, identifies overlaps for various 
    categorical and quantitative biological features of interest. Optionally runs Fishers Exact Tests and Students t-tests.
    
    Parameters:
    - skip_tests (bool): If True, skips the Fisher's Exact Test calculations.
    '''
    if workers > 1:
        pandarallel.initialize(nb_workers=workers, progress_bar=False)
        apply_func = lambda df, func: df.parallel_apply(func)
    else:
        apply_func = lambda df, func: df.progress_apply(func)
    
    CRISPR_compensation, corum_with_symbols, ebi2, humap2, regulons, biogrid_direct, all_screened_pairs, scores_filtered, all_pairs = dfs_for_cat_overlap
    string_physical, biogrid, conservation_scores, prot_hl, interface_df, gtex_data = dfs_for_quant_overlap

    print('Getting datasets ready to calculate overlaps...')
    all_unique_pairs = all_pairs.drop_duplicates(subset="sorted_gene_pair").set_index("sorted_gene_pair")
    bronze_standard_SLs = set(CRISPR_compensation[CRISPR_compensation.SL == True].sorted_gene_pair)
    bronze_standard_tested_pairs = set(CRISPR_compensation.sorted_gene_pair)
    strict_comb_SLs = set(all_screened_pairs[all_screened_pairs.strict_comb_hit].sorted_gene_pair)
    strict_comb_pairs_tested = set(all_screened_pairs[all_screened_pairs.n_screens_tested > 1].sorted_gene_pair)
    lenient_comb_SLs = set(all_screened_pairs[all_screened_pairs.lenient_comb_hit].sorted_gene_pair)
    lenient_comb_pairs_tested = set(all_screened_pairs.sorted_gene_pair)

    proteins_in_corum = set(corum_with_symbols.gene_symbol)
    proteins_in_ebi_complexportal = set(ebi2.gene_symbol)
    proteins_in_humap2 = set(humap2.gene_names)
    overlap_df_output = deepcopy(results)
    overlap_df_output['sorted_gene_pair'] = overlap_df_output['gene_pair'].apply(lambda x: '_'.join(sorted(x.split('_'))))
    essgenes = set(identify_essential_genes(scores_filtered))
    complexes_with_essential_subunits = [
        complex for complex in corum_with_symbols.complex.unique() 
        if essgenes.intersection(corum_with_symbols[corum_with_symbols.complex == complex].gene_symbol)
    ]
    proteins_in_essential_complexes = set(corum_with_symbols[corum_with_symbols.complex.isin(complexes_with_essential_subunits)].gene_symbol)
    
    def check_coregulated(paralog_pair):
        gene1, gene2 = paralog_pair.split('_')
        gene1_tfs = regulons[regulons.target == gene1].tf.to_list()
        gene2_tfs = regulons[regulons.target == gene2].tf.to_list()
        if len(np.intersect1d(gene1_tfs, gene2_tfs)) > 0:
            return True
        else:
            return False
        
    print("Calculating overlap with coreg...")
    regulons['tf_gene_pair'] = regulons['tf'] + '_' + regulons['target']
    overlap_df_output['coreg_tf'] = overlap_df_output['gene_pair'].apply(lambda x: check_coregulated(x))
    
    print("Calculating overlap with direct_interaction...")
    biogrid_direct['gene_pair'] = biogrid_direct['gene_A'] + '_'+  biogrid_direct['gene_B']
    biogrid_direct['sorted_gene_pair'] = biogrid_direct['gene_pair'].apply(lambda x: '_'.join(sorted(x.split('_'))))
    interacting_pairs = list(biogrid_direct.sorted_gene_pair.unique())
    overlap_df_output['direct_interaction'] = overlap_df_output['sorted_gene_pair'].apply(lambda x: x in interacting_pairs)

    print("Calculating overlap with bronze_standard_SL...")
    overlap_df_output["bronze_standard_SL"] = apply_func(
        overlap_df_output["sorted_gene_pair"],
        lambda x: x in bronze_standard_SLs if x in bronze_standard_tested_pairs else np.nan
    )
    print("Calculating overlap with strict_comb_hit...")
    overlap_df_output["strict_comb_hit"] = apply_func(
        overlap_df_output["sorted_gene_pair"],
        lambda x: x in strict_comb_SLs if x in strict_comb_pairs_tested else np.nan
    )
    print("Calculating overlap with lenient_comb_hit...")
    overlap_df_output["lenient_comb_hit"] = apply_func(
        overlap_df_output["sorted_gene_pair"],
        lambda x: x in lenient_comb_SLs if x in lenient_comb_pairs_tested else np.nan
    )
    print("Calculating overlap with in_PC_CORUM...")
    overlap_df_output["in_PC_CORUM"] = apply_func(
        overlap_df_output["sorted_gene_pair"],
        lambda x: any(gene in proteins_in_corum for gene in x.split("_"))
    )
    print("Calculating overlap with in_PC_CORUM_essential...")
    overlap_df_output["in_PC_CORUM_essential"] = apply_func(
        overlap_df_output["sorted_gene_pair"],
        lambda x: any(gene in proteins_in_essential_complexes for gene in x.split("_"))
    )
    print("Calculating overlap with in_PC_CORUM_both...")
    overlap_df_output["in_PC_CORUM_both"] = apply_func(
        overlap_df_output["sorted_gene_pair"],
        lambda x: check_if_both_proteins_in_same_complex(x, corum_with_symbols)
    )
    print("Calculating overlap with in_PC_EBI...")
    overlap_df_output["in_PC_EBI"] = apply_func(
        overlap_df_output["sorted_gene_pair"],
        lambda x: any(gene in proteins_in_ebi_complexportal for gene in x.split("_"))
    )
    print("Calculating overlap with in_PC_humap...")
    overlap_df_output["in_PC_humap"] = apply_func(
        overlap_df_output["sorted_gene_pair"],
        lambda x: any(gene in proteins_in_humap2 for gene in x.split("_"))
    )
    
    print("Mapping closest_pair...")
    all_unique_pairs.reset_index(inplace = True)
    overlap_df_output = overlap_df_output.merge(all_unique_pairs[["sorted_gene_pair", "closest"]], on = 'sorted_gene_pair', how = 'left').rename(
        columns = {'closest':'closest_pair'})
    
    print("Mapping famsize2...")
    overlap_df_output = overlap_df_output.merge(all_unique_pairs[["sorted_gene_pair", "famsize2"]], on = 'sorted_gene_pair', how = 'left')
    all_unique_pairs.set_index('sorted_gene_pair', inplace  = True)
    
    all_tested = set(results.gene_pair)
    all_comp_results = set(results[results.compensation].gene_pair)
    all_cl_results = set(results[results.collateral_loss].gene_pair)
    
    print("Calculating overlap with other_dir_comp...")
    overlap_df_output["other_dir_comp"] = apply_func(
        overlap_df_output["gene_pair"],
        lambda x: f"{x.split('_')[1]}_{x.split('_')[0]}" in all_comp_results 
                  if f"{x.split('_')[1]}_{x.split('_')[0]}" in all_tested else np.nan
    )
    
    print("Calculating overlap with other_dir_cl...")
    overlap_df_output["other_dir_cl"] = apply_func(
        overlap_df_output["gene_pair"],
        lambda x: f"{x.split('_')[1]}_{x.split('_')[0]}" in all_cl_results 
                  if f"{x.split('_')[1]}_{x.split('_')[0]}" in all_tested else np.nan
    )
    
    print("Removing duplicate gene pairs...")
    overlap_df_output = overlap_df_output.drop_duplicates(subset='gene_pair')
    
    output_dfs = [overlap_df_output]

    if not skip_tests:
        print("Preparing interaction datasets for FET analysis...")
        interaction_datasets_list = [
            'bronze_standard_SL', 'strict_comb_hit', 'lenient_comb_hit', 'in_PC_CORUM',
            'in_PC_CORUM_essential', 'in_PC_CORUM_both', 'in_PC_EBI', 'in_PC_humap',
            'closest_pair', 'famsize2', 'other_dir_comp', 'other_dir_cl', 'direct_interaction',
            'coreg_tf'
        ]
        interaction_datasets_list = [x for x in interaction_datasets_list if x in overlap_df_output.columns]
        FETs_df_output_gene_pair = pd.DataFrame()
        FETs_df_output_gene_pair_sorted = pd.DataFrame()
        
        for hit_type in ["compensation", "collateral_loss"]:
            for interaction_dataset in interaction_datasets_list:
                print(f"Running FET for {hit_type} and {interaction_dataset}...")
                FETs_df_output_gene_pair = pd.concat([
                    FETs_df_output_gene_pair, 
                    run_FET(interaction_dataset, hit_type, 'gene_pair', overlap_df_output)
                ])
                FETs_df_output_gene_pair_sorted = pd.concat([
                    FETs_df_output_gene_pair_sorted, 
                    run_FET(interaction_dataset, hit_type, 'sorted_gene_pair', overlap_df_output)
                ])
                
        print("Processing FET results...")
        FETs_df_output_processed_genepair = process_FET(FETs_df_output_gene_pair)
        FETs_df_output_processed_sortedgenepair = process_FET(FETs_df_output_gene_pair_sorted)
        output_dfs += [FETs_df_output_processed_genepair, FETs_df_output_processed_sortedgenepair]

    if runwith == 'prot':
        print("Constructing network from STRING...")
        stringnet_physical = build_ppi_network(string_physical)
        
        print("Constructing network from Biogrid...")
        biogridnet = build_ppi_network(biogrid)
        
        print("Calculating STRING degree centrality...")
        centrality_string_physical = nx.degree_centrality(stringnet_physical)
        
        print("Calculating Biogrid degree centrality...")
        centrality_bg = nx.degree_centrality(biogridnet)

        print("Preparing annotated results dataframe...")
        annotated_results = deepcopy(results)
        
        def add_up(df):
            # Hardcoded choices based on previous selections
            hardcoded_choices = {
                'CDKN2A': 'P42771',
                'PRKCH': 'P24723',
                'SCRIB': 'Q14160',
                'MACF1': 'Q9UPN3',
                'RAB34': 'Q9BZG1',
                'RABGAP1L': 'Q5R372',
                'SLC35A4': 'Q96G79',
                'TOR2A': 'Q5JU69'
            }
            genemap = pd.read_csv('../data/general/geneidmap_sep24.txt', sep = '\t').rename(columns= {'Approved symbol':'gene_name', 'UniProt ID(supplied by UniProt)':'upid'})[['gene_name', 'upid']]
            genemap = genemap.dropna().set_index('gene_name')
            mapdict = genemap.to_dict()['upid']

            def map_protein(protein):
                if protein in hardcoded_choices:
                    return hardcoded_choices[protein]
                uniprot_ids = mapdict.get(protein, "")
                if "," in uniprot_ids:
                    return uniprot_ids.split(",")[0].strip()
                return uniprot_ids
            
            print("Mapping proteins to UniProt IDs...")
            df['A1_uniprot'] = df['A1'].apply(map_protein)
            df['A2_uniprot'] = df['A2'].apply(map_protein)
            return df

        def get_lengths(uniprot_ids):
            results = []
            total = len(uniprot_ids)
            
            for i, uid in enumerate(tqdm(uniprot_ids, desc="Fetching protein lengths")):
                try:
                    r = requests.get(f"https://rest.uniprot.org/uniprotkb/{uid}", params={"format": "json"})
                    if r.status_code == 200:
                        data = r.json()
                        results.append({
                            "uniprot_id": uid,
                            "protein_length": data.get("sequence", {}).get("length")
                        })
                    else:
                        print(f"Failed to fetch {uid}: HTTP {r.status_code}")
                    
                    if i % 10 == 0 and not hasattr(tqdm, "_instances"):
                        print(f"Progress: {i+1}/{total} ({((i+1)/total)*100:.1f}%)")
                        
                    time.sleep(0.3)
                except Exception as e:
                    print(f"Error: {uid} - {e}")
            
            return pd.DataFrame(results)

        print("Adding UniProt IDs to annotated results...")
        annotated_results = add_up(annotated_results)
        all_upids = list(np.union1d(annotated_results['A1_uniprot'].to_list(), annotated_results['A2_uniprot'].to_list()))

        print("Checking for protein length data...")
        if True:
        #if not os.path.exists('../data/for_overlap/prot_lengths.csv'):
            print("Fetching protein lengths from UniProt...")
            lengths = get_lengths(all_upids)
            lengths.to_csv('../data/for_overlap/prot_lengths.csv') # Uncomment below to save time: but UniProt IDs are different for HAP1 and CPTAC pairs so this is default to avoid errors
        #else:
        #    print("Loading protein lengths from file...")
        #    lengths = pd.read_csv('../data/for_overlap/prot_lengths.csv')
            
        print("Mapping protein lengths...")
        lendict = lengths.set_index('uniprot_id').to_dict()['protein_length']
        annotated_results['A1_length'] = annotated_results['A1_uniprot'].apply(lambda x: lendict[x])
        annotated_results['A2_length'] = annotated_results['A2_uniprot'].apply(lambda x: lendict[x])

        print("Mapping interface residues...")
        interface_dict = interface_df.set_index('protein').to_dict()['interface_residues']
        annotated_results['A1_interface'] = annotated_results['A1'].apply(lambda x: interface_dict[x] if x in interface_dict.keys() else np.nan)
        annotated_results['A2_interface'] = annotated_results['A2'].apply(lambda x: interface_dict[x] if x in interface_dict.keys() else np.nan)
        
        print("Calculating interface percentages...")
        annotated_results['A1_perc_at_int'] = (annotated_results['A1_interface'] / annotated_results['A1_length']) * 100
        annotated_results['A2_perc_at_int'] = (annotated_results['A2_interface'] / annotated_results['A2_length']) * 100

        print("Mapping protein half-lives...")
        hl_dict = prot_hl.set_index('gene_name').to_dict()['hl_logged']

        def annotate_hl(df):
            print("Extracting genes from gene pairs...")
            df['A1'] = df['gene_pair'].apply(lambda x: x.split('_')[0])
            df['A2'] = df['gene_pair'].apply(lambda x: x.split('_')[1])
            
            print("Mapping half-life values...")
            df['A1_hl'] = df['A1'].apply(lambda x: hl_dict[x] if x in hl_dict.keys() else np.nan)
            df['A2_hl'] = df['A2'].apply(lambda x: hl_dict[x] if x in hl_dict.keys() else np.nan)
            return df
        
        annotated_results = annotate_hl(annotated_results)
        
        print("Mapping protein abundance data...")
        gtex_dict = gtex_data.set_index('gene_name').to_dict()['mean']

        def annotate_with_abundance(df):
            def annotate_gene_pair(gene_pair):
                A1, A2 = gene_pair.split('_')
                A1_abundance = gtex_dict.get(A1, np.nan)
                A2_abundance = gtex_dict.get(A2, np.nan)
                return pd.Series({
                    'A1_abundance': A1_abundance,
                    'A2_abundance': A2_abundance
                })
            
            print("Adding protein abundance data...")
            result = df['gene_pair'].apply(annotate_gene_pair)
            return pd.concat([df, result], axis=1)
        
        annotated_results = annotate_with_abundance(annotated_results)

        print("Creating sorted gene pairs for network analysis...")
        annotated_results['sorted_gene_pair'] = annotated_results['gene_pair'].apply(lambda x: '_'.join(sorted(x.split('_'))))
        
        print("Calculating STRING physical Jaccard index...")
        annotated_results['string_jaccard_physical'] = apply_func(
            annotated_results['sorted_gene_pair'],
            lambda x: jaccard_index(stringnet_physical, x)
        )
        
        print("Calculating Biogrid Jaccard index...")
        annotated_results['bg_jaccard'] = apply_func(
            annotated_results['sorted_gene_pair'],
            lambda x: jaccard_index(biogridnet, x)
        )
        
        print("Calculating conservation score...")
        annotated_results['conservation_score'] = apply_func(
            annotated_results['A2'],
            lambda x: conservation_scores.loc[x, 'total'] if x in conservation_scores.index else np.nan
        )
        
        print("Calculating median lost gene neighbour essentiality...")
        annotated_results['A2_median_neighbour_ess'] = apply_func(
            annotated_results['A2'],
            lambda x: get_ppi_ess(x, stringnet_physical, scores_filtered)
        )
        
        print("Calculating STRING physical degree centrality...")
        annotated_results['string_dc_physical'] = apply_func(
            annotated_results['A2'],
            lambda x: centrality_string_physical[x] if x in centrality_string_physical else np.nan
        )
        
        print("Calculating Biogrid degree centrality...")
        annotated_results['biogrid_dc'] = apply_func(
            annotated_results['A2'],
            lambda x: centrality_bg[x] if x in centrality_bg else np.nan
        )
        
        print("Calculating sequence identity...")
        annotated_results['sequence_identity'] = apply_func(
            annotated_results['sorted_gene_pair'],
            lambda x: all_unique_pairs.loc[x, 'min_seq_id'] if x in all_unique_pairs.index else np.nan
        )
        
        print("Calculating family size...")
        annotated_results['family_size'] = apply_func(
            annotated_results['sorted_gene_pair'],
            lambda x: all_unique_pairs.loc[x, 'family_size'] if x in all_unique_pairs.index else np.nan
        )
        
        print("Determining hit types...")
        annotated_results['category'] = annotated_results.apply(
            lambda x: get_hit_type(x['compensation'], x['collateral_loss']), 
            axis=1
        )
        
        if not skip_tests:
            print("Preparing variable names for t-tests...")
            namedict = {
                'A2_median_neighbour_ess': 'Median lost gene neighbour essentiality',
                'bg_jaccard': 'Biogrid Jaccard index',
                'biogrid_dc': 'Biogrid degree centrality',
                'conservation_score': 'Conservation score',
                'family_size': 'Family size',
                'sequence_identity': 'Sequence identity',
                'string_dc_physical': 'STRING physical degree centrality',
                'string_jaccard_physical': 'STRING physical Jaccard index',
                'coreg_tf':'Coregulated by same transcription factor',
                'A1_length':'Compensating protein length',
                'A2_length':'Lost protein length', 
                'A1_perc_at_int':'% compensating protein residues\never at interface',
                'A2_perc_at_int':'% lost protein residues\never at interface', 
                'A1_hl':'Log10 compensating\nprotein half life',
                'A2_hl':'Log10 lost protein\nhalf life',
                'A1_dis':'% compensating protein\ndisordered', 
                'A2_dis':'% lost protein\ndisordered',
                'A1_abundance':'Mean compensating\nprotein abundance\nin GTEX', 
                'A2_abundance':'Mean lost protein\nabundance in GTEX'
            }

            print("Running multiple t-tests with unique gene pairs...")
            results_ttests = multiple_ttests_with_correction(
                annotated_results,
                ['string_jaccard_physical', 'bg_jaccard', 'conservation_score', 
                 'sequence_identity', 'family_size'],  drop = True
            )
            
            print("Running multiple t-tests without dropping duplicates...")
            results_ttests_nodrop = multiple_ttests_with_correction(
                annotated_results,
                ['A1_length', 'A2_length', 'A1_perc_at_int', 'A2_perc_at_int',
                 'A2_median_neighbour_ess', 'biogrid_dc', 'string_dc_physical',
                 'A1_hl', 'A2_hl', 'A1_abundance', 
                 'A2_abundance'],  drop = False # Direction is important here so we dont consider only unique gene pairs
            )
            
            print("Combining t-test results...")
            results_ttests_full = pd.concat([results_ttests, results_ttests_nodrop])
            results_ttests_full['Variable'] = results_ttests_full['colname'].apply(lambda x: namedict[x])
            output_dfs.append(results_ttests_full)

        print("Calculating other direction compensation overlap...")
        annotated_results["other_dir_comp"] = apply_func(
            annotated_results["gene_pair"],
            lambda x: f"{x.split('_')[1]}_{x.split('_')[0]}" in all_comp_results 
                  if f"{x.split('_')[1]}_{x.split('_')[0]}" in all_tested else np.nan
        )

        print("Calculating other direction collateral loss overlap...")
        annotated_results["other_dir_cl"] = apply_func(
            annotated_results["gene_pair"],
            lambda x: f"{x.split('_')[1]}_{x.split('_')[0]}" in all_cl_results 
                    if f"{x.split('_')[1]}_{x.split('_')[0]}" in all_tested else np.nan
        )

        print("Adding annotated results to output dataframes...")
        output_dfs.append(annotated_results)
    
    print("Overlap calculation completed.")
    return output_dfs
