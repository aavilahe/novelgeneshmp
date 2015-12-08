#Functions for cleaning subject variables

## Superkeys for HMP tables

#' DCM superkey
DCM_sk = c('PROT',
           'PROTSEG',
           'SITE',
           'RANDSID',
           'MEDNUM'
           )

#' DTPDHXDVD superkey
DTPDHXDVD_sk = c('PROT',
                 'PROTSEG',
                 'SITE',
                 'RANDSID',
                 'VISNO'
                 )

#' DSU superkey
DSU_sk = c('PROT',
           'PROTSEG',
           'SITE',
           'RANDSID'
           )

#' DEM superkey
DEM_sk = c('PROT',
           'protseg',  # Note lowercase!
           'SITE',
           'RANDSID'
           )

#' Extracts medication information from \code{DCM_df}
#'
#' Extract medication information from DCM data.frame
#' per patient (\code{RANDSID}) per visit (\code{VISNO})
#'
#' @param DCM_df A data.frame with \code{DCM} contents
#' @param DTPDHXDVD A data.frame with \code{DTP DHX DVD} contents
#' @return A data.frame with cleaned DCM data.
#'          Rows are samples indexed by a superkey
#'          Columns are various DCMCODEs or the superkey
#' @seealso \code{\link{DCM_sk}}
clean_dcm = function(DCM_df, DTPDHXDVD_df){
    # For each visit figure out if medication was being taken

    DCM_df = left_join(DTPDHXDVD_df %>%
                           select(one_of(DTPDHXDVD_sk), study_day),
                       DCM_df %>%
                           select(one_of(DCM_sk), DCMCODE, MSTART, MSTOP)
                       ) %>%
                filter(MSTART <= study_day & (study_day <= MSTOP | is.na(MSTOP))) %>%
                select(-study_day, -MSTART, -MSTOP)

    # Join again to get DCMCODE = NA (and MEDNUM = NA)
    # for patients that were not on any medication for a particular visit
    DCM_df = DCM_df %>% full_join(DTPDHXDVD_df %>%
                                   select(one_of(DTPDHXDVD_sk))
                                  )

    # Turn DCMCODE into a bunch of binary variables for each DCMCODE value
    DCM_df = DCM_df %>% mutate(hasDCM = TRUE,
                               DCMCODE = paste('DCMCODE_', DCMCODE, sep = '')
                               ) %>%
                        spread(DCMCODE, hasDCM, fill = FALSE)

    # Summarize medication presence/absence by patient per visit (collapse MEDNUM)
    DCM_df = DCM_df %>% group_by_(.dots = c('VISNO',
                                            grep('MEDNUM', DCM_sk,
                                                 invert = TRUE, value = TRUE
                                                 )
                                            )
                                  ) %>%
                        summarise_each(funs(any),
                                       starts_with('DCMCODE_')
                                       ) %>%
                        ungroup()

    # Enforce types
    DCM_df %<>% mutate(PROT = as.integer(PROT),
                       PROTSEG = as.character(PROTSEG),
                       SITE = as.character(SITE),
                       RANDSID = as.integer(RANDSID),
                       VISNO = as.character(VISNO)
                       ) %>%
                mutate_each(funs(as.logical), starts_with('DCMCODE_'))
    return(DCM_df)
}

#' Extracts diet information from \code{DSU_df}
#'
#' Extract diet information from DSU data.frame
#' per patient (\code{RANDSID}).
#'
#' @param DSU_df A data.frame with \code{DSU} contents
#' @return A data.frame with cleaned DSU data.
#'          Rows are samples indexed by a superkey
#'          Columns are DSUBFED, DSUDIET, or the superkey
#' @seealso \code{\link{DSU_sk}}
clean_dsu = function(DSU_df){
    # We only want DSUBFED and DSUDIET

    DSU_df = DSU_df %>% select(one_of(DSU_sk), DSUBFED, DSUDIET) %>%
                        mutate(PROT = as.integer(PROT),           # Enforce types
                               PROTSEG = as.character(PROTSEG),
                               SITE = as.character(SITE),
                               RANDSID = as.integer(RANDSID)
                               ) %>%
                        mutate_each(funs(as.character), DSUBFED, DSUDIET)
    return(DSU_df)
}


#' Extracts demographic information from \code{DEM_df}
#'
#' Extract demographic information from DEM data.frame
#' per patient (\code{RANDSID}). Currently only checks
#' if patient was born in USA/Canada or not.
#'
#' @param DEM_df A data.frame with \code{DEM} contents
#' @return A data.frame with cleaned DEM data.
#'          Rows are samples indexed by a superkey
#'          Columns are MURICA or the superkey
#' @seealso \code{\link{DEM_sk}}
clean_dem = function(DEM_df){
    # We only want BRTHCTRY, and we're going to change it to USA/Canada or not
    # Don't trust data dictionary's "coded values"; double check with BRTHCTRY_C

    DEM_df = DEM_df %>% mutate(MURICA = BRTHCTRY == 63) %>%
                        select(one_of(DEM_sk), MURICA) %>%
                        mutate(PROT = as.integer(PROT),           # Enforce types
                               protseg = as.character(protseg),
                               SITE = as.character(SITE),
                               RANDSID = as.integer(RANDSID)
                               )
    return(DEM_df)
}

#' Extracts lifestyle information from \code{DTPDHXDVD_df}
#'
#' Extract lifestyle information from DTPDHXDVD data.frame
#' per patient (\code{RANDSID}) per (\code{VISNO}).
#' Currently checks for Bachelors degree, student status, smoker status,
#' and BMI.
#'
#' @param DTPDHXDVD_df A data.frame with \code{DTPDHXDVD} contents
#' @return A data.frame with cleaned DTPDHXDVD data.
#'          Rows are samples indexed by a superkey
#'          Columns are lifestyle variables or the superkey
#' @seealso \code{\link{DTPDHXDVD_sk}}
clean_dtpdhxdvd = function(DTPDHXDVD_df){
    # Gets Education, Occupation
    # Splits BMI into three categoris
    # Gets binary variable for cigarette smoker status

    DTPDHXDVD_df = DTPDHXDVD_df %>% select(one_of(DTPDHXDVD_sk),
                                           DVDEDLVL, DVDOCPTN,
                                           DTPBMI,
                                           DVDTOBC, DVDTOBSP
                                           )
    # Education is split into above BS or not
    # Occupation is split into student or not
    # BMI is split into lean, overweight, obese
    # Smoker is a tobacco user that uses cigarettes
    DTPDHXDVD_df = DTPDHXDVD_df %>% mutate(EDLVL_BS = ifelse(is.na(DVDEDLVL), NA,
                                                             DVDEDLVL %in% 8:11
                                                             ),
                                           OCPTN_ST = DVDOCPTN == 42,
                                           BMI_CAT = ifelse(DTPBMI < 25,
                                                             "lean",
                                                             ifelse(DTPBMI < 30,
                                                                    "overweight",
                                                                    "obese"
                                                                    )
                                                             ),
                                           SMOKER = DVDTOBC == 1 &
                                                        (DVDTOBSP == 1 | DVDTOBSP == 3)
                                           ) %>%
                                    select(one_of(DTPDHXDVD_sk),
                                           EDLVL_BS, OCPTN_ST, BMI_CAT, SMOKER
                                           )
    # Enforce types
    DTPDHXDVD_df = DTPDHXDVD_df %>%
                        mutate(PROT = as.integer(PROT),
                               PROTSEG = as.character(PROTSEG),
                               SITE = as.character(SITE),
                               RANDSID = as.integer(RANDSID),
                               VISNO = as.character(VISNO)
                               )
    return(DTPDHXDVD_df)
}

#' Extracts RANDSID, VISNO, SN from GTV
#'
#' Strips leading zeros from primary visits \code{VISNO = (1, 2, 3)}
#'
#' @param GTV_df A data.frame with GTV contents
#' @return A cleaned GTV data.frame with RANDSID, VISNO, SN mapping
clean_gtv = function(GTV_df){
    GTV_df = GTV_df %>% select(SN, RANDSID, VISNO) %>% unique() %>%
                        # Strip leading 0 for integer-like VISNOs (01 -> 1; 01S -> 01S)
                        mutate(VISNO = sub('^0(\\d+)$', '\\1', as.character(VISNO)),
                               RANDSID = as.integer(RANDSID),  # Enforce types
                               SN = as.integer(SN)
                               )
    return(GTV_df)
}
