# Functions for loading and formatting raw FUnkSFAM abundance data files

#' Convert read counts to binary presence-absence.
#'
#' @param cnts_df A data.frame loaded by \code{\link{load_counts}}
#' @param threshold An integer. Counts must be greater than this for presence.
#' @return A presence-absence data.frame.
#'          Rows are samples indexed by SRSID (first column)
#'          Columns except the first are FUnkSFAMs
counts_to_presence = function(cnts_df, threshold = 4){
    FUNKIDs = grep('SRSID', colnames(cnts_df), invert = TRUE, value = TRUE)
    pres_df = cnts_df %>% mutate_each_(funs(ifelse(. > threshold, TRUE, FALSE)),
                                 FUNKIDs
                                 )
    return(pres_df)
}

#' Load counts and return presence-absence data.frame.
#' @param config_l A list parsed from config file
#' @return A presence-absence data.frame.
#' @seealso \code{\link{counts_to_presence}}
prepare_abundance = function(config_l){
    cnts_df = load_counts(config_l$counts_fn)
    pres_df = counts_to_presence(cnts_df, config_l$count_threshold)
    return(pres_df)
}

#' Make mapping from (RANDSID, VISNO) to SRSID including body (sub)site information.
#' @param config_l A list parsed from config file
#' @return A data.frame with columns:
#'          SN, RANDSID, VISNO, SRSID, HMP_BodySite, HMP_BodySubsite
prepare_map = function(config_l){
    GTV_df = clean_gtv(load_pheno(config_l$GTV_fn))
    s2v = load_srs2visno(config_l$srs2visno_fn)
    s2r = load_srs2randsid(config_l$srs2randsid_fn)
    bodysites_df = load_bodysites(config_l$project_catalog_fn)
    srs_map = bind_rows(inner_join(GTV_df, s2v), filter(s2r, SRSID %ni% s2v$SRSID)) %>% distinct(SRSID)
    srs_map %<>% left_join(bodysites_df)
    return(srs_map)
}
