import pandas as pd


def get_pairs_to_test(pairs, runwith, cndf, only_run_HAP1_pairs=False,
                      seq_id_thresh=0.2, fam_size_thresh=20):
    """
    Filter down from Ensembl's full set of paralog pairs to the set that we'll
    eventually test, while keeping track of how many pairs are dropped at each
    step.
    Drop genes without prot/trans/resid and CN data, apply sequence identity
    threshold, family size threshold, drop paralog pairs on same chromosome
    Drop paralog pairs where either member of the pair is on a sex chromosome 
    """
    print('Getting paralog pairs to test...')
    common_genes = list(cndf.gene_name.unique())  # Already filtered, saves adding another variable to params
    pairs_in_quant_and_CN = pairs[pairs['gene_pair'].apply(
        lambda x: check_if_both_genes_in_set(x, common_genes))]
    pairs_in_quant_and_CN_minseq = pairs_in_quant_and_CN[
        pairs_in_quant_and_CN.min_seq_id >= seq_id_thresh]
    pairs_in_quant_and_CN_minseq_maxfam = pairs_in_quant_and_CN_minseq[
        pairs_in_quant_and_CN_minseq.family_size <= fam_size_thresh]
    pairs_in_quant_and_CN_minseq_maxfam_diffchr = (
        pairs_in_quant_and_CN_minseq_maxfam[
            pairs_in_quant_and_CN_minseq_maxfam.same_chr == False
        ]
    )
    pairs_to_return = pairs_in_quant_and_CN_minseq_maxfam_diffchr[
        pairs_in_quant_and_CN_minseq_maxfam_diffchr.sex_chr == False]

    filtering_dict = {
        'All pairs': len(pairs),
        f'Pairs with {runwith} and copy num data': len(pairs_in_quant_and_CN),
        f'Min seq ID >= {seq_id_thresh}': len(pairs_in_quant_and_CN_minseq),
        f'Max fam size <= {fam_size_thresh}': len(pairs_in_quant_and_CN_minseq_maxfam),
        f'A1 and A2 on different chromosomes': len(pairs_in_quant_and_CN_minseq_maxfam_diffchr),
        f'No sex chromosome genes': len(pairs_to_return)
    }
    
    if only_run_HAP1_pairs:
        HAP1_results = pd.read_csv('../output/HAP1/HAP1_paralogtest_results.csv')
        pairs_to_return = pairs_to_return[
            pairs_to_return.gene_pair.isin(HAP1_results.gene_pair.to_list())
        ]
        filtering_dict['Tested in HAP1 analysis'] = len(pairs_to_return)
    return filtering_dict, pairs_to_return


def process_pairs_to_test(pairs_to_test, lost, filtering_dict, n_lost_limit=20):
    """
    Process pairs to test by dropping genes that are lost too few times.
    """
    n_A2_lost = lost.sum(axis=1).reset_index().rename(
        columns={'index': 'A2', 'gene_name': 'A2', 0: 'n_A2_lost'}
    )
    df = pd.merge(pairs_to_test, n_A2_lost, how='inner')
    pairs_to_test_processed = df[(df.n_A2_lost >= n_lost_limit)]
    filtering_dict[f'Pairs where A2 is lost at least {n_lost_limit} times'] = len(pairs_to_test_processed)
    return pairs_to_test_processed, filtering_dict


def filter_unbacked_losses(pairs_to_test_processed, results_A2, filtering_dict):
    """
    Process pairs to test by dropping pairs where loss of a paralog is not
    associated with a significant drop in its own protein abundance.
    """
    backed_losses = set(results_A2[results_A2.backed].A2.unique())
    pairs_to_test_processed_backed = pairs_to_test_processed[
        pairs_to_test_processed.A2.isin(backed_losses)
    ]
    filtering_dict['Tested A2s'] = len(results_A2)
    filtering_dict[f'A2s whose loss results in a drop in A2 protein'] = len(backed_losses)
    return pairs_to_test_processed_backed, filtering_dict

def check_if_both_genes_in_set(gene_pair, gene_set):
    """
    For a given gene pair and a given gene set, check if both genes in the pair are in the set.
    """
    A1, A2 = gene_pair.split('_')
    return (A1 in gene_set) and (A2 in gene_set)
