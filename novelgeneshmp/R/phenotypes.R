# Functions for loading and cleaning subject metadata

#' Load all phenotype data in config file
#'
#' @param config_l
#' @return A list of cleaned data.frames
#' @export
load_and_clean_phenotypes = function(config_l){
    # Load
    DCM_raw       = load_pheno(config_l$DCM_fn)
    DSU_raw       = load_pheno(config_l$DSU_fn)
    DEM_raw       = load_pheno(config_l$DEM_fn)
    DTPDHXDVD_raw = load_pheno(config_l$DTPDHXDVD_fn)

    # Clean
    phenos_l = list(DCM_df        = clean_dcm(DCM_raw, DTPDHXDVD_raw),
                    DSU_df        = clean_dsu(DSU_raw),
                    DEM_df        = clean_dem(DEM_raw),
                    DTPDHXDVD_df  = clean_dtpdhxdvd(DTPDHXDVD_raw)
                    )
    return(phenos_l)
}

#' Get phenotype variable names from column names
#'
#' NOTE: Names and name regexes are hardcoded
#'
#' @param column_names All variable names that may include some phenotype names
#' @return a vector of phenotype variable names
#' @export
get_phenos_from_colnames = function(column_names){
    #column_names = colnames(df)
    pheno_names = column_names[column_names %in%
                               c('MURICA',
                                 'BMI_CAT',
                                 'SMOKER',
                                 'OCPTN_ST',
                                 'EDLVL_BS',
                                 'DSUDIET',
                                 'DSUBFED',
                                 grep('^DCMCODE_', column_names, value = TRUE)
                                 )
                               ]
    return(pheno_names)
}

#' Get VISNO values that can be mapped to samples (SRSIDs).
#'
#' @param srs_map A mapping of SRSID to RANDSID, VISNO and other metadata and
#'                identifiers
#' @return A vector of mappable VISNO values
#' @export
get_mappable_visnos = function(srs_map){
    mappable_visnos = unique(na.omit(srs_map$VISNO))  # get mappable VISNOs (e.g. 1, 2, 3)
    cat('DBG: Returning mappable VISNOs\n', file = stderr())
    cat(paste(mappable_visnos, '\n'), file = stderr())
    return(mappable_visnos)
}

#' Inspect phenotype metadata mapping.
#'
#' Phenotype data can be mapped to samples of a particular VISNO or to all VISNOs.
#' However only abundance data was available for VISNO = 1,2,3 and not the
#' supplementary visits.
#'
#' Prints what percentage or counts of samples remain after mapping each way
#' for the user to inspect.
#'
#' @param phenos_l A named list of data.frames with phenotype data
#' @param mappable_visnos A vector of VISNO values that can be mapped to
#' @export
check_mappable_phenos = function(phenos_l, mappable_visnos){
    cat('DBG: ----------------------------------------------------\n', file = stderr())
    cat('DBG: Manual check for mapping reads to phenotype by VISNO\n', file = stderr())
    cat('DBG: ----------------------------------------------------\n', file = stderr())

    cat('DBG: DCM: Check if there are more VISNO-mappable than -unmappable samples\n\n', file = stderr())
    sink(stderr())
    print(
          phenos_l$DCM_df %>%
            mutate(IS_MAPPABLE = VISNO %in% mappable_visnos) %>%
            group_by(IS_MAPPABLE) %>%
            tally()
        )
    cat('\n', file = stderr())
    sink()
    cat('DBG: DTPDHXDVD: Check % !is.na of VISNO-mappable vs -unmappable samples\n', file = stderr())
    cat('DBG: DTPDHXDVD: for OCPTN_ST, EDLVL_BS, BMI_CAT, and SMOKER\n\n', file = stderr())
    sink(stderr())
    print(
          phenos_l$DTPDHXDVD_df %>%
            mutate(IS_MAPPABLE = VISNO %in% mappable_visnos) %>%
            group_by(IS_MAPPABLE) %>%
            summarise_each(funs({sum(!is.na(.))/length(.)}), OCPTN_ST, EDLVL_BS, BMI_CAT, SMOKER)
          )
    cat('\n', file = stderr())
    sink()
}

#' Joins VISNO-mappable phenotypes from \code{check_mappable_phenos}.
#'
#' NOTE: Hardcoded behavior after inspecting printed output of
#'       \code{\link{check_mappable_phenos}}
#'
#' @param phenos_l A list of named phenotype data.frames
#' @param mappable_visnos A vector of mappable VISNO values
#' @return A data.frame of VISNO-mappable phenotypes
#'         Columns: RANDSID, VISNO, SITE, [PHENONAME1, ...]
#' @seealso \code{\link{check_mappable_phenos}}
#' @export
get_visno_mappable_pheno_df = function(phenos_l, mappable_visnos){
    visno_mappable_phenos = c(grep('^DCMCODE_', colnames(phenos_l$DCM_df), value = TRUE),
                             'OCPTN_ST', 'EDLVL_BS', 'SMOKER'
                             )

    cat('DBG: visno-mappable phenos: ', file = stderr())
    cat(paste(visno_mappable_phenos, '\n'), file = stderr())

    list_of_dfs = list(phenos_l$DCM_df, phenos_l$DTPDHXDVD_df)
    visno_pheno_df = Reduce(full_join, list_of_dfs) %>%
                        select(one_of(visno_mappable_phenos),
                               RANDSID, VISNO, SITE
                               ) %>%
                        filter(VISNO %in% mappable_visnos)
    return(visno_pheno_df)
}

#' Joins VISNO-unmappable phenotypes from \code{check_mappable_phenos}.
#'
#' NOTE: Hardcoded behavior after inspecting printed output of
#'       \code{\link{check_mappable_phenos}}
#'
#' Unmappable phenotypes are mapped to samples from ALL VISNOs.
#'
#' @param phenos_l A list of named phenotype data.frames
#' @param mappable_visnos A vector of mappable VISNO values
#' @return A data.frame of VISNO-unmappable phenotypes
#'         Columns: RANDSID, SITE, [PHENONAME1, ...]
#' @seealso \code{\link{check_mappable_phenos}}
#' @export
get_visno_unmappable_pheno_df = function(phenos_l, mappable_visnos){
    # Joins list_of_dfs, selects (RANDSID, visno_unmappable_phenos)

    visno_unmappable_phenos = c('MURICA', 'DSUBFED', 'DSUDIET',
                                'BMI_CAT'
                                )

    cat('DBG: visno-unmappable phenos: ', file = stderr())
    cat(paste(visno_unmappable_phenos, '\n'), file = stderr())

    list_of_dfs = list(phenos_l$DEM_df, phenos_l$DSU_df,
                       phenos_l$DTPDHXDVD_df %>% filter(VISNO %ni% mappable_visnos) %>%
                           select(RANDSID, BMI_CAT)
                       )

    nonvisno_pheno_df = Reduce(full_join, list_of_dfs) %>%       # join multiple dfs
                        select(one_of(visno_unmappable_phenos),  # select only relevant columns
                               RANDSID, SITE
                               ) %>%
                        group_by(RANDSID) %>%                    # drop rows with NAs in duplicate RANDSIDs
                            do({if(nrow(.) > 1){na.omit(.)} else {.}}) %>%
                            ungroup()
    return(nonvisno_pheno_df)
}
