pick_tests = function(column_names){
    ff_names = grep('^X', column_names, value = TRUE)
    #not_ph_names = unique(c(ff_names,
    #                        'HMP_BodySite',
    #                        'HMP_BodySubsite',
    #                        'SRSID',
    #                        'SN',
    #                        DSU_sk,
    #                        DEM_sk,
    #                        DCM_sk,
    #                        DTPDHXDVD_sk
    #                        )
    #                      )
    #ph_names = column_names[column_names %ni% not_ph_names]
    ph_names = get_phenos_from_colnames(column_names)
    ff_ph = setNames(expand.grid(ff_names, ph_names),
                     c('FUNKID', 'PHENONAME')
                     ) %>%
            mutate_each(funs(as.character))

    #warn("DBG: SUBSAMPLING FOR TESTING PURPOSES\n")
    #if(nrow(ff_ph) > 20){
    #    ff_ph %<>% sample_n(20)
    #}

    return(ff_ph)
}

get_res_df = function(df, ff_name, ph_name){
    the_formula = as.formula(paste(ff_name, '~ SITE +', ph_name))
    the_model = try(
                    glm(formula = the_formula, data = df,
                        family = 'binomial'
                        ),
                    silent = TRUE  # suppresses error messages
                    )
    if('try-error' %ni% class(the_model)){
        summary_df = model_to_summary_df(the_model)
        stats_df = get_model_stats(the_model, ff_name, ph_name)
        res_df = summary_df %>% rowwise() %>% do({cbind(., stats_df)})
        return(res_df)
    } else {
        warn(paste(":::> There was an error testing:",
                   ff_name, ph_name, "\n",
                   collapse = " "
                   )
             )
        return(data.frame())
    }
}

glm_loop = function(df){
    df %<>% prefilter_by_arbitrary_statistic(
                calc_PHvariation_per_group(df, group = NULL)
                )

    ff_ph = pick_tests(colnames(df))
    warn(paste0("\nAbout to fit ", nrow(ff_ph), " models!\n"))
    res_df = ff_ph %>% rowwise() %>%
                    do({get_res_df(df, .$FUNKID, .$PHENONAME)})
    return(res_df)
}

do_glm_tests = function(df, group_name){
    res_df = df %>% group_by_(.dots = group_name) %>%
                do({glm_loop(.)}) %>%
                ungroup()
    return(res_df)
}

model_to_summary_df = function(the_model){
    summary_df = summary(the_model)$coefficients %>% as.data.frame() %>%
                add_rownames() %>% tbl_df()
    colnames(summary_df) = c('Coeff_Name', 'Estimate',
                         'Std_Err', 'z_value', 'p_value'
                         )
    return(summary_df)
}

get_model_stats = function(the_model, ff_name, ph_name){
    stats_df = data.frame(PHENONAME = ph_name, FUNKID = ff_name,
                          N_samples = the_model$model %>% nrow(),
                          PH_Entropy_bits = the_model$model[, ph_name] %>% get_entropy(),
                          PH_ArbitraryStatistic = the_model$model[, ph_name] %>% get_arbitrary_statistic(),
                          FF_Entropy_bits = the_model$model[, ff_name] %>% get_entropy(),
                          FF_PercentPresent = the_model$model[, ff_name] %>% get_percent_present()
                          )
    return(stats_df)
}

