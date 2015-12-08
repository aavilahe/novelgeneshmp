# Functions for computing FUnkSFAM and phenotype stats

#' Gets percent presence of FUNKID across samples
#'
#' @param ff_pres A presence-absence logical vector
#' @return A percentage
get_percent_present = function(ff_pres){
    return(sum(ff_pres) / length(ff_pres))
}

#' Gets entropy (James-Stein shrinkage estimator) for a categorical column
#'
#' Requires \code{entropy} package
#'
#' @param cat_column A categorical vector
#' @return Entropy in bits
get_entropy = function(cat_column){
    return(entropy::entropy(table(cat_column),
                            verbose = FALSE,
                            method = 'shrink',
                            unit = 'log2'
                            )
           )
}

#' Gets ArbitraryStatistic for a categorical column
#'
#' I don't know what to call this.
#'
#' Checks that there is more than 1 value with more than 4 counts
#'
#' @param cat_column A categorical vector
#' @return ArbitraryStatistic
get_arbitrary_statistic = function(cat_column){
    return(sum(table(cat_column) > 4) > 1)
}

#' Get number of samples per group (body site or subsite)
#'
#' @param df A data.frame in `tidy` format
#' @param group Either 'HMP_BodySite' or 'HMP_BodySubsite'
#' @return A data.frame
calc_nsamples_per_group = function(df, group){
    nsamp_df = df %>% group_by_(.dots = group) %>% tally()
    return(nsamp_df)
}

#' Get FUnkSFAM variation per group
#'
#' @inheritParams calc_nsamples_per_group
#' @return A data.frame wih Entropy_bits and PercentPresent per FUNKID per group
calc_FFvariation_per_group = function(df, group){
    variation_df = df %>% gather(FUNKID, Present, starts_with('X')) %>%
                        group_by_(.dots = c(group, 'FUNKID')) %>%
                        summarise_each(funs(Entropy_bits = get_entropy,
                                            PercentPresent = get_percent_present
                                            ),
                                       Present
                                       ) %>%
                        ungroup()
    return(variation_df)
}

#' Get phenotype variation per group
#'
#' @inheritParams calc_nsamples_per_group
#' @return A data.frame with Entropy_bits and ArbitraryStatistic per PHENONAME
#'         per group
calc_PHvariation_per_group = function(df, group){
    column_names = colnames(df)
    #pheno_names = column_names[column_names %in%
    #                           c('MURICA',
    #                             'BMI_CAT',
    #                             'SMOKER',
    #                             'OCPTN_ST',
    #                             'EDLVL_BS',
    #                             'DSUDIET',
    #                             'DSUBFED',
    #                             grep('^DCMCODE_', column_names, value = TRUE)
    #                             )
    #                           ]
    pheno_names = get_phenos_from_colnames(column_names)
    variation_df = df %>% gather(PHENONAME, PHENOVAL, one_of(pheno_names)) %>%
                        group_by_(.dots = c(group, 'PHENONAME')) %>%
                        summarise_each(funs(ArbitraryStatistic = get_arbitrary_statistic,
                                            Entropy_bits = get_entropy
                                            ),
                                       PHENOVAL
                                       ) %>%
                        ungroup()

    return(variation_df)
}

#' DEPRECATED. Use \code{\link{calc_PHvariation_per_group}}
calc_PHas_per_group = function(df, group){
    column_names = colnames(df)
    #pheno_names = column_names[column_names %in%
    #                           c('MURICA',
    #                             'BMI_CAT',
    #                             'SMOKER',
    #                             'OCPTN_ST',
    #                             'EDLVL_BS',
    #                             'DSUDIET',
    #                             'DSUBFED',
    #                             grep('^DCMCODE_', column_names, value = TRUE)
    #                             )
    #                           ]
    pheno_names = get_phenos_from_colnames(column_names)
    PHas_df = df %>% group_by_(.dots = group) %>%
                        summarise_each(funs(get_arbitrary_statistic), one_of(pheno_names)) %>%
                        gather(PHENONAME, ArbitraryStatistic, one_of(pheno_names))
    return(PHas_df)
}

