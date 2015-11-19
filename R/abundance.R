counts_to_presence = function(cnts_df, threshold = 4){
    # Convert "Counts" table to binary "Presence/absence"
    # Present if Counts > threshold (default 4); absence otherwise

    FUNKIDs = grep('SRSID', colnames(cnts_df), invert = TRUE, value = TRUE)
    pres_df = cnts_df %>% mutate_each_(funs(ifelse(. > threshold, TRUE, FALSE)),
                                 FUNKIDs
                                 )
    return(pres_df)
}

prepare_abundance = function(config_l){
    # Load
    cnts_df = load_counts(config_l$counts_fn)
    pres_df = counts_to_presence(cnts_df, config_l$count_threshold)
    return(pres_df)
}

prepare_map = function(config_l){
    # Map to RANDSID, VISNO
    GTV_df = clean_gtv(load_pheno(config_l$GTV_fn))
    s2v = load_srs2visno(config_l$srs2visno_fn)
    s2r = load_srs2randsid(config_l$srs2randsid_fn)
    bodysites_df = load_bodysites(config_l$project_catalog_fn)
    srs_map = bind_rows(inner_join(GTV_df, s2v), filter(s2r, SRSID %ni% s2v$SRSID)) %>% distinct(SRSID)
    srs_map %<>% left_join(bodysites_df)
    return(srs_map)
}
