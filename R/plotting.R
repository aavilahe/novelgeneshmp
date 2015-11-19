
get_quantiles_df = function(column){
    quantile_df = data.frame(PercentileValue = quantile(column)) %>%
                        add_rownames(var = "Percentile")
    return(quantile_df)
}

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

