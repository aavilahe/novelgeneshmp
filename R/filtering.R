# Filters

#' Removes tests where FUnkSFAM abundance is mostly constant
#'
#' Currently hardcoded to keep if presence-absence entropy
#' is greater than 75%-ile for a particular group (e.g. HMP_BodySite,
#' HMP_Bodysubsite)
#'
#' @param res_df A results data.frame computed by \code{\link{do_glm_tests()}}
#' @param ff_stats_df A FUnkSFAM variation statistics data.frame computed by
#'                    \code{\link{calc_FFvariation_per_group()}}
#' @param group Either 'HMP_Bodysite' or 'HMP_Bodysubsite'
#' @return A filtered res_df
filter_results_by_FFentropy_per_group = function(res_df, ff_stats_df, group){
    warn("DBG: Hardcoded filter: FFentropy > 75%-ile\n")
    q_df = ff_stats_df %>% group_by_(.dots = group) %>%
            do({get_quantiles_df(.$Entropy_bits)}) %>%
            filter(Percentile == '75%')
    keep_df = ff_stats_df %>% left_join(q_df) %>%
                    filter(Entropy_bits > PercentileValue) %>%
                    select(one_of(group, 'FUNKID'))
    filtered_res_df = keep_df %>% inner_join(res_df)
    return(filtered_res_df)
}

#' Removes tests where FUnkSFAM annotation is high or phylogeny shallow
#'
#' `Span` must be 6 or 7 and \code{N_annotations} must be less than 10
#'
#' @param res_df A results data.frame computed by \code{\link{do_glm_tests()}}
#' @param config_l Configuration list read by \code{\link{load_config}}
#' @seealso load_config
#' @return A filtered res_df
filter_results_by_funksfam_annotation = function(res_df, config_l){
    # Loads span_6_7_annotations.tsv, keeps families with less than 10 annotations

    res_df %<>% separate(FUNKID, c("Xdigit", "NoX_FFID"),
                         sep = '_', remove = FALSE
                         ) %>%
                mutate_each(funs(as.character), NoX_FFID, FUNKID)

    span67_fn = config_l$span67_fn
    keep_df = load_tsv(span67_fn) %>%
                    select(FUNKSFAM, N_annotations) %>%
                    mutate(NoX_FFID = as.character(FUNKSFAM)) %>%
                    select(-FUNKSFAM) %>%
                        # FUNKSFAMs should be unique, keep highest number of annotations
                        group_by(NoX_FFID) %>%
                        summarise(N_annotations = max(N_annotations)) %>%
                        ungroup() %>%
                    filter(N_annotations < 10)
    filtered_res_df = keep_df %>% inner_join(res_df) %>%
                        select(-Xdigit, -NoX_FFID)
    return(filtered_res_df)
}

#' Removes tests where the subject variable being tested doesn't vary enough
#'
#' "ArbitraryStatistic" for subject variable (aka medical metadata) is computed
#' by \code{\link{calc_PHvariation_per_group()}}
#'
#' @param res_df A results data.frame computed by \code{\link{do_glm_tests()}}
#' @return A filtered res_df
filter_results_by_per_test_arbitrary_statistic = function(res_df){
    # This used to be prefiltered to decrease compute time and glm errors
    filtered_res_df = res_df %>% filter(PH_ArbitraryStatistic)
    return(filtered_res_df)
}

#' Removes DCMCODEs not considered for the study
#'
#' Keeps only "interesting" DCMCODEs (see Methods in accompanying manuscript)
#' - DCMCODE_NA means no medication was taken for given sample
#' - DCMCODEs of interest (all DCMCODEs were tested):
#'     - 03 (antacids)
#'     - 08 (antibiotics)
#'     - 13 (contraceptives)
#'     - 16 (GI meds)
#'     - 18 (hormones/steroids)
#'
#' @param res_df A results data.frame computed by \code{\link{do_glm_tests()}}
#' @return A filtered res_df
filter_results_by_interesting_DCMCODE = function(res_df){

    warn("DBG: Hardcoded filter: DCMCODE_XX\n")
    interesting_DCMCODEs = paste0('DCMCODE_', c(
                                      8,
                                      3,
                                      13,
                                      16,
                                      18
                                      )
                                  )
    filtered_res_df = res_df %>% filter(!grepl('DCMCODE_', PHENONAME) |
                                        PHENONAME %in% interesting_DCMCODEs
                                        )
    return(filtered_res_df)
}

#' Removes tests where Coeff_Name doesn't start with PHENONAME
#'
#' Removes unneeded tests (e.g. SITE, (Intercept))
#'
#' @param res_df A results data.frame computed by \code{\link{do_glm_tests()}}
#' @return A filtered res_df
filter_unneeded_tests = function(res_df){
    # Remove unneeded tests
    #
    # e.g. (Intercept), SITE92WAU

    res_df = res_df %>% rowwise() %>%
                filter(grepl(paste0('^', PHENONAME), Coeff_Name)) %>%
                ungroup()  # (p.adjust fails without this)
    return(res_df)
}

#' Runs all post-testing filters
#'
#' See \code{vignette('example-run')}
#'
#' @param res_df A results data.frame computed by \code{\link{do_glm_tests()}}
#' @param ff_stats_df A FUnkSFAM variation statistics data.frame computed by
#'                    \code{\link{calc_FFvariation_per_group()}}
#' @param group Either 'HMP_Bodysite' or 'HMP_Bodysubsite'
#' @return A filtered res_df
filter_results_all = function(res_df, ff_stats_df, group, config_l){
    res_df %<>% filter_unneeded_tests() %>%
                filter_results_by_FFentropy_per_group(ff_stats_df, group) %>%
                filter_results_by_funksfam_annotation(config_l) %>%
                filter_results_by_per_test_arbitrary_statistic() %>%
                filter_results_by_interesting_DCMCODE()
    return(res_df)
}

#' Prefilter sub-data.frames prior to model fitting
#'
#' Similar to \code{\link{filter_results_by_per_test_arbitrary_statistic}}
#'
#' @param df A data.frame to fit a series of glm models
#' @param ph_stats_df Precomputed "ArbitraryStatistic"s for each sample
#'                      (does not account for samples that may be dropped if
#'                       corresponding FUnkSFAM presence is NA)
#' @return A filtered df
prefilter_by_arbitrary_statistic = function(df, ph_stats_df){
    warn("\nPrefiltering by ArbitraryStatistic\n")
    ph_names = get_phenos_from_colnames(colnames(df))
    ph_keep = ph_stats_df %>%
                group_by(PHENONAME) %>%
                summarise(ArbitraryStatistic = any(ArbitraryStatistic)) %>%
                filter(ArbitraryStatistic) %>%
                extract2("PHENONAME") %>%  # `[[`
                as.character()  # or else select() doesn't always work
    prefiltered_df = df %>% select(-one_of(ph_names),
                                   one_of(ph_keep[ph_keep %in% ph_names])
                                   )
    return(prefiltered_df)
}
