import numpy as np
import pandas as pd 
import statsmodels.formula.api as smf
from pandarallel import pandarallel
from tqdm import tqdm

def run_test(pair, lost, quants, lin, how, cn_df, runwith):
    """
    Fits a linear model to explain the abundance of one protein with the loss status of a gene 
    (either its own gene i.e. A2-A2 tests or another i.e paralogous gene, A1-A2). 
    """
    indlist = ['gene_pair', 'A1', 'A2', 'ols_p', 'ols_coef', 'n_A2_lost', 'lost_mean_quant', 'other_mean_quant']
    
    # Always include gene pair information
    result = pd.Series([pair.gene_pair, pair.A1, pair.A2] + [np.nan] * 5, index=indlist)
    if (pair.A1 in quants.columns.to_list()) and (pair.A2 in lost.index.to_list()):
        quants_sliced = quants.loc[:, pair[how]].reset_index().rename(columns={'index':'sample_ID', 'Sample':'sample_ID','':'sample_ID', pair[how]:'quant'})
        A2_info = lost.loc[pair.A2,:].reset_index().rename(columns={'index':'sample_ID', '':'sample_ID', pair.A2:'A2_lost'})
        df = pd.merge(quants_sliced, A2_info).dropna(subset=['quant', 'A2_lost'])
        df = df.merge(lin, on = 'sample_ID', how = 'inner') # Always merge with lineage information as this can always be a confounder
        if how == 'A1': # We're running A1-A2 tests for paralog abundance
            cn_sliced = cn_df[cn_df.gene_name == pair['A1']].drop(columns = 'gene_name')
            df = df.merge(cn_sliced, on = 'sample_ID', how = 'inner')
        if df.shape[0]!=0:
            if how == 'A2': # We're running self i.e. A2-A2 tests to confirm that loss is associated with a drop in protein abundance
                ols_results = smf.ols('quant ~ C(A2_lost) + C(lineage)', data=df).fit() # Cancer type/lineage is a covariate for the self-test
            elif how == 'A1':
                if runwith in ['prot', 'trans']:
                    model = 'quant ~ C(A2_lost) + C(lineage) + CNV' # Cancer type/lineage and self-copy number are covariates for the paralog test
                elif runwith == 'prot_residual':
                    model = 'quant ~ C(A2_lost)'# Those have both been regressed out in the residuals so we don't need them here
                ols_results = smf.ols(model, data=df).fit()
            ols_pval = ols_results.pvalues['C(A2_lost)[T.True]'] # t-statistic of the model is extracted as a p-value
            ols_coef = ols_results.params['C(A2_lost)[T.True]'] # OLS coefficient of the model is extracted
            if (pair.gene_pair == 'FAM219B_FAM219A'):
                 print(quants_sliced.head())
                 print(A2_info.head())
                 print(df.head())
                

            if not np.isnan(ols_pval): # We successfully ran the test
                result = pd.Series([pair.gene_pair, pair.A1, pair.A2, ols_pval, ols_coef, df[df.A2_lost == True].shape[0], 
                                    df[df.A2_lost == True].quant.mean(), df[df.A2_lost == False].quant.mean()],
                                    index = indlist)
    return result             

def run_tests_for_df(how, workers, pairs_to_test_df, lost, quants, lin, cn_df, runwith):
      """
      Run the linear regression test for each gene pair in pairs_to_test as per requirements.
      """
      print(f'Running {how}-A2 regressions. This step can take a while...')
      if how == 'A2':
          pairs_to_test_df = pairs_to_test_df.drop_duplicates(subset = 'A2') # Only run one test per lost protein
      if workers > 1: # If parallel processing has been enabled-- by default it is not, pandarallel used if it is
            pandarallel.initialize(nb_workers=workers, progress_bar=True) # Set to False to avoid screen clutter, can turn on to monitor progress 
            res_raw = pairs_to_test_df.parallel_apply(run_test, args = (lost, quants, lin, how, cn_df, runwith), axis = 1)
      else:
            res_raw = pairs_to_test_df.progress_apply(run_test, args = (lost, quants, lin, how, cn_df, runwith), axis = 1)
      if how == 'A2':
           res_raw = res_raw.drop(columns = ['gene_pair', 'A1'])
      print(f'Done with {how}-A2 tests!')
      return res_raw
