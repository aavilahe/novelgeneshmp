# Function to load and save various files and file formats

#' Load a csv file with a header
#'
#' @param fn A filename
#' @return A data.frame
#' @export
load_csv = function(fn){
    df = tbl_df(
                read.csv(fn, header = TRUE, sep = ',',
                    na.string = c('', NA)
                    )
                )
    return(df)
}

#' Load a tab-separated value file with a header
#'
#' @param fn A filename
#' @return A data.frame
#' @export
load_tsv = function(fn){
    df = tbl_df(
                read.table(fn, header = TRUE, sep = '\t',
                    na.string = c('', NA)
                    )
                )
    return(df)
}


#' Load raw phenotype data (DSU, DCM, etc...)
#'
#' Checks extension if \code{\link{load_tsv}} or
#' \code{\link{load_csv}} should be called
#'
#' @param fn A filename
#' @return A data.frame with the phenotype data
#' @export
load_pheno = function(fn){
    if(grepl('.txt$', fn)){
        # Load txt file
        pheno_df = load_tsv(fn)
    } else {
        # Load csv file
        pheno_df = load_csv(fn)
    }
    return(pheno_df)
}

#' Load abundance data (as counts)
#'
#' Loads and transposes counts csv
#' where each row is a FUnkSFAM. The first column is a FUNKSFAM ID,
#' other columns are counts for each SRSID
#'
#' @param fn A filename
#' @return A data.frame with each row a sample and FUnkSFAMs as columns
#' @export
load_counts = function(fn){
    cnts_df = load_csv(fn)
    colnames(cnts_df)[1] = 'FUNKID'
    cnts_df$X.1 = NULL  # There was an extra column to delete
    cnts_df %<>% gather(key = 'SRSID', value = 'Counts', -FUNKID) %>%
              select(FUNKID, SRSID, Counts) %>%
              mutate(FUNKID = paste('X', FUNKID, sep = '')) %>%  # "X" preprended to colnames
              spread(FUNKID, Counts) %>%
              mutate(SRSID = as.character(SRSID)) %>%    # Enforce types
              mutate_each(funs(as.integer), starts_with('X'))  # Enforce types
    return(cnts_df)
}


#' Load various specific maps
#'
#' Reads HMP project catalog in csv format
#' and selects Sequence.Read.Archive.ID, HMP.Isolation.Body.Site, HMP.Isolation.Body.Subsite
#' renames columns to SRSID, HMP_BodySite, HMP_BodySubsite
#'
#' @param fn A filename
#' @return A data.frame with SRSIDs mapped to body sites and subsites
#' @export
load_bodysites = function(fn){
    srs_bodysites_df = load_csv(fn) %>%
                            rename(SRSID = Sequence.Read.Archive.ID,
                                   HMP_BodySite = HMP.Isolation.Body.Site,
                                   HMP_BodySubsite = HMP.Isolation.Body.Subsite
                                   ) %>%
                            select(SRSID, HMP_BodySite, HMP_BodySubsite) %>%
                            mutate_each(funs(as.character))  # Enforce types
    return(srs_bodysites_df)
}

#' Loads a mapping found in ppAll_V35_map.txt
#'
#' @param fn A filename
#' @return A data.frame mapping of SRSID to VISNO (and SN for good measure)
#' @export
load_srs2visno = function(fn){
    srs2visno_df = load_tsv(fn) %>%
                        select(SRSID = SRS_SampleID,
                               VISNO = VisitNo,
                               SN
                               ) %>%
                        unique() %>%
                        mutate(SRSID = as.character(SRSID),  # Enforce types
                               VISNO = as.character(VISNO),
                               SN = as.integer(SN)
                               )
    return(srs2visno_df)
}


#' Loads a mapping from SRSID to RANDSID
#'
#' Maps metagenomic sample IDs to randomized patient IDs
#'
#' @param fn A filename
#' @return A data.frame mapping RANDSID to SRSID
#' @export
load_srs2randsid = function(fn){
    srs2randsid_df = load_tsv(fn) %>%
                        rename(RANDSID = dbGaP.SubjID) %>%  # Fixes column name
                        mutate(SRSID = as.character(SRSID),  # Enforce types
                               RANDSID = as.integer(RANDSID)
                               )
    return(srs2randsid_df)
}

#' Saves a data.frame to a tab-separated value file
#'
#' Writes column names as header, No rownames, No quotes.
#'
#' @param df A data.frame to save
#' @param fn A filename to save the df to
#' @seealso load_tsv, load_csv
#' @export
save_tsv = function(df, fn){
    write.table(df, file = fn,
                sep = '\t',
                col.names = TRUE, row.names = FALSE,
                quote = FALSE
                )
}

#' Formats final results
#'
#' Selects interesting results to print, sorts by adjusted and
#' unadjusted p-values
#'
#' @param res_df A data.frame of results
#' @return A formatted res_df
#' @export
format_final_results = function(res_df){
    # choose columns and order to report in table
    res_df %<>% select(matches('^HMP_BodyS(ite|ubsite)$'),
                       FUNKID, PHENONAME, Coeff_Name, p_adj_fdr, p_value,
                       N_samples,
                       FF_Entropy_bits, FF_PercentPresent,
                       PH_Entropy_bits,
                       Estimate, Std_Err, z_value
                       ) %>%
                arrange(p_adj_fdr, p_value)
    return(res_df)
}
