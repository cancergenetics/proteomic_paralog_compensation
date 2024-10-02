import pandas as pd
import numpy as np

def load_and_process_data(runwith='prot', overlap_only=False):
    """
    Load and preprocess the required datasets.
    
    Parameters:
    - runwith: Specifies which dataset to load ('prot', 'trans', or 'prot_residual')
    - overlap_only: If True, only load and return the overlap dataframes
    
    Returns:
    If overlap_only is False:
        ((data_df, cndf, sample_info, all_pairs), dfs_for_cat_overlap, dfs_for_quant_overlap)
    If overlap_only is True:
        (dfs_for_cat_overlap, dfs_for_quant_overlap)
    """
    if not overlap_only:
        print('Loading in data...')
        if runwith == 'prot':
            data_path = '../data/CPTAC/CPTAC_PanCan_proteomics.parquet'
        elif runwith == 'trans':
            data_path = '../data/CPTAC/CPTAC_PanCan_transcriptomics.parquet'
        elif runwith == 'prot_residual':
            data_path = '../data/CPTAC/CPTAC_prot_residuals.csv'
        cn_path = '../data/CPTAC/CPTAC_gistic_cn.csv'
        sample_info_path = '../data/CPTAC/CPTAC_sampleinfo.csv'
        
        print('Loading in sample info...')
        sample_info = pd.read_csv(sample_info_path).rename(
            columns={
                "Unnamed: 0": "sample_ID",
                "DepMapID": "sample_ID",
                "Study": "lineage"
            }
        )[["sample_ID", "lineage"]]

        print(f'Loading in {runwith} data...')
        if 'parquet' in data_path:
            data_df = pd.read_parquet(data_path)
        elif 'csv' in data_path:
            data_df = pd.read_csv(data_path, index_col=0, low_memory=False)

        if runwith == 'prot_residual':
            data_df = data_df[['gene_in_sample', 'Residuals']]
            data_df[['gene', 'sample']] = data_df['gene_in_sample'].str.split('_in_', n=1, expand=True)
            data_df = data_df.pivot(index='sample', columns='gene', values='Residuals')
        
        print(f'Data df starts at {data_df.shape}')
        print('Loading in copy number data...')
        cndf = pd.read_csv(cn_path)

        print('Processing datasets...')
        # Keeping genes and samples with copy number data
        data_df.columns.name = ''
        data_df.index.name = 'sample_ID'
        data_df = data_df[data_df.index.isin(sample_info[sample_info.lineage == 'OvCa 2020'].sample_ID.to_list()) == False]
        cndf = cndf.melt(id_vars='gene_name').rename(
            columns={'variable': 'sample_ID', 'value': 'CNV'}
        )
        cndf['gene_in_sample'] = cndf['gene_name'] + '_in_' + cndf['sample_ID']
        cndf.set_index('gene_in_sample', inplace=True)

        common_samples = list(np.intersect1d(cndf.sample_ID.to_list(), data_df.index.to_list()))
        common_genes = list(np.intersect1d(cndf.gene_name.to_list(), data_df.columns.to_list()))
        data_df = data_df[data_df.index.isin(common_samples)][common_genes]
        cndf = cndf[(cndf.sample_ID.isin(common_samples)) & (cndf.gene_name.isin(common_genes))]
        data_df.columns = data_df.columns.astype(str)
        data_df = data_df.dropna(axis=1, thresh=(len(data_df) // 2))
        common_samples = list(np.intersect1d(cndf.sample_ID.to_list(), data_df.index.to_list()))
        common_genes = list(np.intersect1d(cndf.gene_name.to_list(), data_df.columns.to_list()))
        data_df = data_df[data_df.index.isin(common_samples)][common_genes]
        cndf = cndf[(cndf.sample_ID.isin(common_samples)) & (cndf.gene_name.isin(common_genes))]

    print('Loading in files needed to calculate biological overlaps...')
    
    pairs_path = '../data/general/all_pairs_wsexchr.csv'
    all_pairs = pd.read_csv(pairs_path, low_memory=False)
    all_pairs['gene_pair'] = all_pairs['A1'] + '_' + all_pairs['A2']
    all_pairs['sorted_gene_pair'] = all_pairs['gene_pair'].apply(lambda x: '_'.join(sorted(x.split('_'))))
    all_pairs['famsize2'] = (all_pairs['family_size'] == 2)

    def get_gene_name(entid, mapping):
       return mapping.loc[entid, 'gene_symbol'] if entid in mapping.index else np.nan

    CRISPR_compensation = pd.read_csv("../data/for_overlap/1-s2.0-S240547122100329X-mmc5.csv")
    CRISPR_compensation['gene_pair'] = CRISPR_compensation['A1'] + '_' + CRISPR_compensation['A2']
    CRISPR_compensation['sorted_gene_pair'] = CRISPR_compensation['gene_pair'].apply(lambda x: '_'.join(sorted(x.split('_'))))
    corum_with_symbols = pd.read_csv("../data/for_overlap/corum_with_symbols_nov23.csv")
    ebi2 = pd.read_csv("../data/for_overlap/ebi_complex_portal.tsv", sep="\t")
    ebi2["uniprot_ids"] = ebi2["Expanded participant list"].str.split("|")
    ebi2 = ebi2.explode("uniprot_ids")
    ebi2["Uniprot_ID"] = ebi2["uniprot_ids"].str[:-3]
    mapping = pd.read_csv("../data/general/geneidmap_sep24.txt", sep="\t")
    mapping = mapping.dropna().rename(columns={"Approved symbol":"gene_symbol", "UniProt ID(supplied by UniProt)":"Uniprot_ID", 'NCBI Gene ID(supplied by NCBI)':'entrez_id'})
    ebi2 = ebi2.merge(mapping[["gene_symbol","Uniprot_ID"]], on="Uniprot_ID", how="inner")    
    humap = pd.read_csv("../data/for_overlap/humap2_complexes_20200809.csv")
    humap2 = humap[humap.Confidence <= 3].copy()
    humap2["gene_names"] = humap2["genenames"].str.split()
    humap2 = humap2.explode("gene_names")
    all_screened_pairs = pd.read_csv("../data/for_overlap/all_screened_paralog_pairs_25_04_22.csv")
    all_screened_pairs["strict_comb_hit"] = all_screened_pairs.n_screens_SL > 1
    all_screened_pairs["lenient_comb_hit"] = all_screened_pairs.n_screens_SL > 0
    scores_filtered = pd.read_csv('../data/for_overlap/scores_depmap.csv', index_col='cell_line')
    mapping['entrez_id'] = mapping['entrez_id'].astype(int).astype(str)
    mapping = mapping.set_index('entrez_id')[['gene_symbol']]
    scores_filtered.columns = scores_filtered.columns.map(lambda x: get_gene_name(x, mapping)).astype(str)

    dfs_for_cat_overlap = (CRISPR_compensation, corum_with_symbols, ebi2, humap2, all_screened_pairs, scores_filtered, all_pairs)

    # Need the rest of these to calculate quantitative overlaps with t-tests
    biomart_enspmap = pd.read_csv('../data/general/biomart.txt', sep ='\t')
    biomart_enspmap = biomart_enspmap[['Protein stable ID', 'Gene name']]
    biomart_enspmap = biomart_enspmap.rename(columns = {'Protein stable ID':'ENSP_ID', 'Gene name':'gene_name'})
    mapdict = biomart_enspmap.dropna().set_index('ENSP_ID').to_dict()['gene_name']
    
    def map_string(stringdb, map_dict):
          def extract_gene_name(protein_id, map_dict):
              split_id = protein_id.split('9606.')[1] if '9606.' in protein_id else None
              return map_dict.get(split_id, np.nan)
          stringdb['gene_name_A'] = stringdb['protein1'].apply(lambda x: extract_gene_name(x, map_dict))
          stringdb['gene_name_B'] = stringdb['protein2'].apply(lambda x: extract_gene_name(x, map_dict))
          stringdb = stringdb[['gene_name_A', 'gene_name_B']].rename(columns = {'gene_name_A':'gene_A', 'gene_name_B':'gene_B'})
          return stringdb
    
    string_physical = pd.read_csv('../data/for_overlap/9606.protein.physical.links.v12.0.txt', sep = ' ')
    string_physical = map_string(string_physical, map_dict = mapdict)
    biogrid = pd.read_csv('../data/for_overlap/biogrid_full.tsv', sep = '\t')
    biogrid = biogrid[biogrid['Interaction Detection Method'] != 'psi-mi:"MI:0254"(genetic interference)']
    biogrid['gene_A'] = biogrid['Alt IDs Interactor A'].apply(lambda x:x.split('locuslink:')[1].split('|')[0])
    biogrid['gene_B'] = biogrid['Alt IDs Interactor B'].apply(lambda x:x.split('locuslink:')[1].split('|')[0])
    biogrid = biogrid[['gene_A', 'gene_B']]

    conservation_scores = pd.read_csv('../data/for_overlap/ens111_humanGene_phylogeneticProfile_1472Sp.csv')
    conservation_scores = conservation_scores[['gene', 'total']].set_index('gene')
     
    dfs_for_quant_overlap = (string_physical, biogrid, conservation_scores)

    if overlap_only:
        return dfs_for_cat_overlap, dfs_for_quant_overlap
    else:
        return (data_df, cndf, sample_info, all_pairs), dfs_for_cat_overlap, dfs_for_quant_overlap