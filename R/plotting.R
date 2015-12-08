#' Functions to plot FUnkSFAM entropy

#' Gets 0%, 25%, 50%, 75%, 100% quantiles of a column
#'
#' Note quantile and percentile used interchangeably here, oops.
#'
#' @param column A numeric column
#' @return A data.frame with percentiles and percentile values
get_quantiles_df = function(column){
    quantile_df = data.frame(PercentileValue = quantile(column)) %>%
                        add_rownames(var = "Percentile")
    return(quantile_df)
}

#' Makes a histogram of FUnkSFAM entropies per group
#'
#' Uses ggplot2.
#'
#' @param ff_stats_df A data.frame of FUnkSFAM entropies
#'                     calculated by calc_FFvariation_per_group()
#' @param group Either 'HMP_BodySite' or 'HMP_Bodysubsite'
#' @return A ggplot2 object
plot_FFentropy_by_group = function(ff_stats_df, group){
    quantile_df = ff_stats_df %>% group_by_(.dots = group) %>%
                        do({get_quantiles_df(.$Entropy_bits)})

    g = ggplot(ff_stats_df) + geom_histogram(aes(x = Entropy_bits)) +
        geom_vline(data = quantile_df,
                   aes(xintercept = PercentileValue,
                       color = Percentile,
                       linetype = Percentile
                       ),
                   size = 0.5,
                   show_guide = TRUE
                   ) +
        facet_wrap(as.formula(paste0('~', group))) +
        theme_bw() + theme(aspect.ratio = 1) +
        ggtitle("FUnkSFAM Presence-Absence Entropy")
    return(g)
}

