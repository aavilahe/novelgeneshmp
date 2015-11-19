# load file formats
load_csv = function(fn){
    df = tbl_df(
                read.csv(fn, header = TRUE, sep = ',',
                    na.string = c('', NA)
                    )
                )
    return(df)
}

load_tsv = function(fn){
    df = tbl_df(
                read.table(fn, header = TRUE, sep = '\t',
                    na.string = c('', NA)
                    )
                )
    return(df)
}


# load raw phenotype data (DSU, DCM, etc...)
load_pheno = function(fn){
    # Reads csv or txt files based on file name

    if(grepl('.txt$', fn)){
        # Load txt file
        pheno_df = load_tsv(fn)
    } else {
        # Load csv file
        pheno_df = load_csv(fn)
    }
    return(pheno_df)
}

# load abundance data (as counts)
load_counts = function(fn){
    # Loads and transposes counts csv
    #
    # Each row is a funkfam
    # First column is a FUNKSFAM ID, other columns are counts for each SRSID
    #
    # Returns a tbl_df with each row a sample and funkfams as columns

    cnts_df = load_csv(fn)
    colnames(cnts_df)[1] = 'FUNKID'
    cnts_df$X.1 = NULL
    cnts_df %<>% gather(key = 'SRSID', value = 'Counts', -FUNKID) %>%
              select(FUNKID, SRSID, Counts) %>%
              mutate(FUNKID = paste('X', FUNKID, sep = '')) %>%  # "X" preprended to colnames
              spread(FUNKID, Counts) %>%
              mutate(SRSID = as.character(SRSID)) %>%    # Enforce types
              mutate_each(funs(as.integer), starts_with('X'))  # Enforce types
    return(cnts_df)
}


# load various specific maps
load_bodysites = function(fn){
    # Reads HMP project catalog in csv format
    # and selects Sequence.Read.Archive.ID, HMP.Isolation.Body.Site, HMP.Isolation.Body.Subsite
    # renames columns to SRSID, HMP_BodySite, HMP_BodySubsite

    srs_bodysites_df = load_csv(fn) %>%
                            rename(SRSID = Sequence.Read.Archive.ID,
                                   HMP_BodySite = HMP.Isolation.Body.Site,
                                   HMP_BodySubsite = HMP.Isolation.Body.Subsite
                                   ) %>%
                            select(SRSID, HMP_BodySite, HMP_BodySubsite) %>%
                            mutate_each(funs(as.character))  # Enforce types
    return(srs_bodysites_df)
}

load_srs2visno = function(fn){
    # Loads a mapping found in ppAll_V35_map.txt
    #
    # Returns mapping of SRSID to VISNO (and SN for good measure)

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

load_srs2randsid = function(fn){
    # Loads map from SRSID to RANDSID
    # maps metagenomic sample IDs to randomized patient IDs

    srs2randsid_df = load_tsv(fn) %>%
                        rename(RANDSID = dbGaP.SubjID) %>%  # Fixes column name
                        mutate(SRSID = as.character(SRSID),  # Enforce types
                               RANDSID = as.integer(RANDSID)
                               )
    return(srs2randsid_df)
}

save_tsv = function(df, fn){
    write.table(df, file = fn,
                sep = '\t',
                col.names = TRUE, row.names = FALSE,
                quote = FALSE
                )
}

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
