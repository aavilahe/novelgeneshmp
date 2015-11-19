adjust_pvalues = function(res_df){
    # Multiple testing correction

    res_df %<>% mutate(p_adj_fdr = p.adjust(p_value, method = 'fdr'))
    return(res_df)
}
