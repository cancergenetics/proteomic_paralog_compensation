import sys
import pandas as pd
from parse_args import parse_args
from data_preprocessing import load_and_process_data
from filter_paralog_pairs import get_pairs_to_test, process_pairs_to_test, filter_unbacked_losses
from call_gene_loss import call_gene_loss
from regressions_tests import run_tests_for_df
from process_results import process_results
from generate_overlaps import generate_overlaps


def main():
    args = parse_args()
    print(f'\n ******* Running paralog compensation analysis with CPTAC {args.runwith} '
          f'data using {args.nb_workers} workers ******* \n')
    
    processed_data = load_and_process_data(runwith=args.runwith)
    data_df, cndf, sample_info, all_pairs = processed_data[0]
    print(f'We are working with data df at {data_df.shape}')
    filtering_dict, pairs_to_test = get_pairs_to_test(
        pairs=all_pairs,
        runwith=args.runwith,
        cndf=cndf,
        only_run_HAP1_pairs=args.HAP1_overlap
    )
    
    A2_lost = call_gene_loss(
        pairs=pairs_to_test,
        samples=list(cndf.sample_ID.unique()),
        cndf=cndf
    )

    if args.runwith in ['trans', 'prot_residual']:
        try:
            prot_results = pd.read_csv('../output/output_CPTAC/prot/paralog_tests_prot.csv')
            pairs_to_test_processed_backed = prot_results[['gene_pair', 'A1', 'A2', 'n_A2_lost']]
        except FileNotFoundError:
            print("Run the analysis with proteomic data first as only these pairs are tested "
                  "with transcriptomic and residual data.")
            sys.exit(1)
    else:
        pairs_to_test_processed, filtering_dict = process_pairs_to_test(
            pairs_to_test=pairs_to_test,
            lost=A2_lost,
            filtering_dict=filtering_dict
        )
        res_raw_A2 = run_tests_for_df(
            how='A2',
            workers=int(args.nb_workers),
            pairs_to_test_df=pairs_to_test_processed,
            lost=A2_lost,
            quants=data_df,
            lin=sample_info,
            cn_df=cndf,
            runwith=args.runwith
        )
        results_A2, filtering_dict = process_results(
            how='A2',
            res_raw=res_raw_A2,
            pairs=pairs_to_test_processed,
            filtering_dict=filtering_dict
        )
        pairs_to_test_processed_backed, filtering_dict = filter_unbacked_losses(
            pairs_to_test_processed=pairs_to_test_processed,
            results_A2=results_A2,
            filtering_dict=filtering_dict
        )
    print(f'Starting with {len(pairs_to_test_processed_backed)} pairs to test')
    res_raw_A1A2 = run_tests_for_df(
        how='A1',
        workers=int(args.nb_workers),
        pairs_to_test_df=pairs_to_test_processed_backed,
        lost=A2_lost,
        quants=data_df,
        lin=sample_info,
        cn_df=cndf,
        runwith=args.runwith
    )
    results_A1A2, filtering_dict = process_results(
        how='A1A2',
        res_raw=res_raw_A1A2,
        pairs=pairs_to_test_processed_backed,
        filtering_dict=filtering_dict
    )
    
    filtering_df = pd.DataFrame(list(filtering_dict.items()), columns=['Step', 'Pairs remaining'])

    dfs_for_cat_overlap = processed_data[1]
    dfs_for_quant_overlap = processed_data[2]

    output_dfs = generate_overlaps(
        dfs_for_cat_overlap,
        dfs_for_quant_overlap,
        results_A1A2,
        args.runwith,
        args.nb_workers
    )

    if args.HAP1_overlap == False:
        if args.runwith == 'prot':
            overlap_df_output, FETs_genepair, FETs_sortedgenepair, quant_ttest_overlaps, annotated_res = output_dfs
            print(f'In total, we identify {len(results_A2[results_A2.backed])} proteins out of '
                f'{len(results_A2)} tested that display a drop in protein abundance when hemizygously lost.\n A2-A2 results saved.')
            
            results_A2.to_csv(f'../output/output_CPTAC/{args.runwith}/self_tests_{args.runwith}.csv')
            filtering_df.to_csv(f'../output/output_CPTAC/{args.runwith}/pairs_filtered_at_each_step_{args.runwith}.csv')
            A2_lost.to_csv(f'../output/output_CPTAC/prot/A2_lost.csv')
        else:
            overlap_df_output, FETs_genepair, FETs_sortedgenepair = output_dfs     

    print(f'This allows us to test {len(results_A1A2)} paralog pairs to associate hemizygous '
          f'gene loss with changes in paralog abundance.')
    print(f'We identify {len(results_A1A2[results_A1A2.compensation])} cases of compensation '
          f'and {len(results_A1A2[results_A1A2.collateral_loss])} cases of collateral loss.')
    print('Unless limiting to HAP1 pairs- biological overlaps have been identified with '
          'Fishers Exact Tests and t-tests, see output dir')
    print('Saving results and exiting...')

    output_filename = f'../output/output_CPTAC/{args.runwith}/paralog_tests_{args.runwith}'
    if args.HAP1_overlap:
        output_filename += '_onlyHAP1pairs'
    output_filename += '.csv'

    results_A1A2.to_csv(output_filename)

    if args.HAP1_overlap == False:
        overlap_df_output.to_csv(f'../output/output_CPTAC/{args.runwith}/categorical_overlaps_{args.runwith}.csv')
        #FETs_genepair.to_csv(f'../output/output_CPTAC/{args.runwith}_/categorical_FETs_directionalgenepairs_{args.runwith}.csv') 
        FETs_sortedgenepair.to_csv(f'../output/output_CPTAC/{args.runwith}/categorical_FETs_uniquegenepairs_{args.runwith}.csv')

        if args.runwith == 'prot':
            quant_ttest_overlaps.to_csv(f'../output/output_CPTAC/{args.runwith}/quantitative_ttests_foroverlaps_prot.csv')
            annotated_res.to_csv(f'../output/output_CPTAC/{args.runwith}/all_quantoverlaps_prot.csv')
    print('*****')


if __name__ == "__main__":
    main()