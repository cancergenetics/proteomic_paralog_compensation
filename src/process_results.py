import pandas as pd
from statsmodels.stats.multitest import fdrcorrection


def process_results(how, res_raw, pairs, filtering_dict, FDRlimit=0.05,
                    pvallimit=0.05, n_lost_limit=20):
    FDRlimit = 0.05
    pvallimit = 0.05
    # We can end up with some tests run with fewer than the minimum number of
    # samples due to NAs, drop these
    res_processed = res_raw[res_raw.n_A2_lost >= n_lost_limit]

    if how == 'A1A2':
        res_processed = res_processed.dropna(subset=["ols_p"]).merge(
            pairs[['gene_pair', 'A1', 'A2']]
        )
        # Apply Benjamini Hochberg multiple testing correction
        res_processed = res_processed.assign(
            p_adj=fdrcorrection(res_processed.ols_p.values)[1]
        )
        res_processed = res_processed.assign(
            compensation=res_processed.apply(
                lambda x: ((x.p_adj < FDRlimit) & (x.ols_p < pvallimit) &
                           (x.ols_coef > 0)),
                axis=1
            )
        )
        res_processed = res_processed.assign(
            collateral_loss=res_processed.apply(
                lambda x: ((x.p_adj < FDRlimit) & (x.ols_p < pvallimit) &
                           (x.ols_coef < 0)),
                axis=1
            )
        )
        # Hits are called at an FDR of 5% and uncorrected p-value threshold
        # of 0.05, with OLS coefficient used to determine effect direction

    if how == 'A2':
        res_processed = res_processed.drop_duplicates(subset='A2')
        res_processed = res_processed.assign(
            backed=res_processed.apply(
                lambda x: ((x.ols_p < pvallimit) & (x.ols_coef < 0)),
                axis=1
            )
        )
        # Similarly for the self test using p-value threshold 0.05 by default

    if how == 'A1A2':
        filtering_dict['Paralog pairs tested in the regression'] = len(res_processed)
        filtering_dict[f'Compensation pairs (FDR<{FDRlimit}, pval<{pvallimit})'] = len(
            res_processed[res_processed.compensation]
        )
        filtering_dict[f'Collateral loss pairs (FDR<{FDRlimit}, pval<{pvallimit})'] = len(
            res_processed[res_processed.collateral_loss]
        )
        filtering_dict['No significant effect'] = len(
            res_processed[(~res_processed.compensation) &
                          (~res_processed.collateral_loss)]
        )
        filtering_dict['Compensation percentage'] = round(
            (len(res_processed[res_processed.compensation]) /
             len(res_processed)) * 100,
            1
        )
        filtering_dict['Collateral loss percentage'] = round(
            (len(res_processed[res_processed.collateral_loss]) /
             len(res_processed)) * 100,
            1
        )

    return res_processed, filtering_dict