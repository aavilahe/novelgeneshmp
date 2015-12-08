#' Multiple testing correction
#'
#' Adds a new column (p_adj_fdr) with fdr adjusted p-values
#'
#' @param res_df A data.frame of results
#' @return res_df with adjust p-values
adjust_pvalues = function(res_df){

    res_df %<>% mutate(p_adj_fdr = p.adjust(p_value, method = 'fdr'))
    return(res_df)
}
