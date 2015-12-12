The `novelgeneshmp` package accompanies supplemental methods for **NOVELGENES\_TITLE\_HERE** by **AUTHORLIST**. It contains functions to facilitate tests for association between subject metadata from the Human Microbiome Project (HMP) and under-annotated protein families investigated in the manuscript (FUnkSFAMs).

Finding associations between FUnkSFAM abundance and subject metadata in Human Microbiome Project samples
--------------------------------------------------------------------------------------------------------

``` r

library(novelgeneshmp)
```

I like using `dplyr` and `magrittr` for manipulating data.frames, and `ggplot2` for plotting.

``` r

library(dplyr)
library(magrittr)
library(ggplot2)
```

### Configuration file

File paths and other options are specified in a `json` formatted configuration file that drives the majority of data 'clean-up' and output.

``` r

# Specify configuration filename
config_fn = system.file("extdata",
                        "novelgenes_configv.json",
                        package = "novelgeneshmp"
                        )

# Uses jsonlite::fromJSON to load the config into a list
config_l = load_config(config_fn)
```

#### Example config:

``` json
{
    "DCM_fn"             : "/path/to/phs000228.v3.pht001187.v3.p1.c1.EMMES_HMP_DCM.HMP.txt",
    "DSU_fn"             : "/path/to/phs000228.v3.pht002156.v1.p1.c1.EMMES_HMP_DSU.HMP.csv",
    "DEM_fn"             : "/path/to/phs000228.v3.pht002158.v1.p1.c1.EMMES_HMP_DEM_ENR.HMP.txt",
    "DTPDHXDVD_fn"       : "/path/to/phs000228.v3.pht002157.v1.p1.c1.EMMES_HMP_DTP_DHX_DVD.HMP.csv",
    "GTV_fn"             : "/path/to/phs000228.v3.pht001193.v3.p1.c1.EMMES_HMP_GTV.HMP.txt",

    "counts_fn"          : "/path/to/HMP_funksfams_counts.csv",

    "project_catalog_fn" : "/path/to/HMP_project_catalog.csv",
    "srs2randsid_fn"     : "/path/to/SRSID2RANDSID.tsv",
    "srs2visno_fn"       : "/path/to/ppAll_V35_map.txt",

    "count_threshold"    : 1,

    "span67_fn"          : "/path/to/span67.tsv",

    "output_dir"         : "./output_dir/",
    "inspect_dir"        : "./inspect_dir/"
}
```

#### Subject metadata

Specify where the subject metadata is. This data is protected and may be requested from [dbGaP](http://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000228.v3.p1). We will not redistribute it.

<table style="width:122%;">
<colgroup>
<col width="26%" />
<col width="22%" />
<col width="73%" />
</colgroup>
<thead>
<tr class="header">
<th align="left">config attribute</th>
<th align="left">data</th>
<th align="left">variables</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left"><code>DCM_fn</code></td>
<td align="left">DCM</td>
<td align="left">Medications</td>
</tr>
<tr class="even">
<td align="left"><code>DSU_fn</code></td>
<td align="left">DSU</td>
<td align="left">Diet</td>
</tr>
<tr class="odd">
<td align="left"><code>DEM_fn</code></td>
<td align="left">DEM</td>
<td align="left">Birthcountry</td>
</tr>
<tr class="even">
<td align="left"><code>DTPDHXDVD_fn</code></td>
<td align="left">DTP DHX DVD</td>
<td align="left">BMI, smoker status, occupation, etc...</td>
</tr>
<tr class="odd">
<td align="left"><code>GTV_fn</code></td>
<td align="left">GTV</td>
<td align="left">Used for <code>RANDSID</code> to (<code>SRSID</code>, <code>VISNO</code>) mapping</td>
</tr>
</tbody>
</table>

#### FUnkSFAM read counts (abundance -&gt; presence-absence)

`counts_fn`:

-   CSV formatted
-   Each row is a FUnkSFAM
-   First column is a FUnkSFAM
-   Other columns are counts for each `SRSID`

`count_threshold`:

-   FUnkSFAMs are considered present if there are more than `count_threshold` read counts

#### Auxiliary data and data mappings

<table style="width:174%;">
<colgroup>
<col width="37%" />
<col width="62%" />
<col width="73%" />
</colgroup>
<thead>
<tr class="header">
<th align="left">config attribute</th>
<th align="left">data</th>
<th align="left">variables</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left"><code>project_catalog_fn</code></td>
<td align="left">HMP project catalog</td>
<td align="left">Body site and subsites</td>
</tr>
<tr class="even">
<td align="left"><code>srs2randsid_fn</code></td>
<td align="left">A mapping from SRSID to RANDSID</td>
<td align="left"><code>SRSID</code>, <code>RANDSID</code></td>
</tr>
<tr class="odd">
<td align="left"><code>srs2visno_fn</code></td>
<td align="left">A mapping from SRSID to VISNO</td>
<td align="left"><code>SRSID</code>, <code>VISNO</code> <em>only contains VISNO = {1, 2, 3}</em></td>
</tr>
<tr class="even">
<td align="left"><code>span67_fn</code></td>
<td align="left">Number of annotations across databases</td>
<td align="left"><code>FUNKID</code>, <code>N_annotations</code></td>
</tr>
</tbody>
</table>

#### Output

These attributes aren't required, but they are handy.

-   `output_dir`: Put results here
-   `inspect_dir`: Put intermediate files here

In these examples, the directory names are expected to end in a `/`. Directories can be created if they don't exist with `dir.create(path, showWarnings = FALSE)`.

### Data preprocessing.

Subject metadata is joined and formatted into two tables. One that contains subject metadata that was collected at each visit, and one that contains subject metadata that is not expected to change across visit (e.g. birth country).

``` r

phenos_l = load_and_clean_phenotypes(config_l)

# Loads counts and convert to presence-abscence
pres_df = prepare_abundance(config_l)
srs_map = prepare_map(config_l)  # maps `SRSID` to (`RANDSID`, `VISNO`)

# get `VISNO`s that are mappable to samples
mappable_visnos = get_mappable_visnos(srs_map)  # Only primary visits in this case

# Inspect printed output for DCMCODE_XX, etc...
check_mappable_phenos(phenos_l, mappable_visnos)

# make phenotype tables, using `check_mappable_phenos` output.
ph_vm = get_visno_mappable_pheno_df(phenos_l, mappable_visnos)     # complicated functions
ph_vu = get_visno_unmappable_pheno_df(phenos_l, mappable_visnos)   # with hardcoded phenotypes


# merge phenotype tables with presence-abscence tables
visno_data_df = Reduce(inner_join, list(pres_df, srs_map, ph_vm))
novisno_data_df = Reduce(inner_join, list(pres_df, srs_map, ph_vu))
```

Because the data is protected, our example proceeds with randomly generated placeholder data.

``` r

# Specify where cleaned data would be
visno_data_fn = system.file("extdata",
                            paste0(config_l$inspect_dir, "visno_mapped_samples_cleaned.tsv"),
                            package = "novelgeneshmp"
                            )
novisno_data_fn = system.file("extdata",
                            paste0(config_l$inspect_dir, "visno_unmapped_samples_cleaned.tsv"),
                            package = "novelgeneshmp"
                            )

# Load the cleaned data
visno_data_df = load_tsv(visno_data_fn)
novisno_data_df = load_tsv(novisno_data_fn)
```

### Variable summaries

Get summaries for phenotypes (subject metadata) and FUnkSFAM abundance (presence-absence).

``` r

# Number of samples per body (sub)site for data that maps to VISNO of origin
# and for data that is not expected to change across VISNOs
nsamples_by_bodysite_df = bind_rows(   # dplyr::bind_rows() is like rbind()
    calc_nsamples_per_group(visno_data_df, 'HMP_BodySite') %>%
        mutate(VISNO_MAPPED = TRUE),
    calc_nsamples_per_group(novisno_data_df, 'HMP_BodySite') %>%
        mutate(VISNO_MAPPED = FALSE)
    )

# FUnkSFAM variation statistics by body (sub)site.
ff_stats_bss_df = calc_FFvariation_per_group(novisno_data_df, 'HMP_BodySubsite')

# Phenotype (subject metadata variable) variation statistics.
ph_stats_bs_df = bind_rows(calc_PHvariation_per_group(visno_data_df, 'HMP_BodySite'),
                           calc_PHvariation_per_group(novisno_data_df, 'HMP_BodySite')
                           )
#> Warning: attributes are not identical across measure variables; they will
#> be dropped
#> Warning in rbind_all(x, .id): Unequal factor levels: coercing to character
```

We can view the summaries:

``` r

head(nsamples_by_bodysite_df)
#> Source: local data frame [6 x 3]
#> 
#>   HMP_BodySite     n VISNO_MAPPED
#>         (fctr) (int)        (lgl)
#> 1           aa    17         TRUE
#> 2           ba    22         TRUE
#> 3           ca    20         TRUE
#> 4           da    21         TRUE
#> 5           ea    22         TRUE
#> 6           aa    85        FALSE
head(ph_stats_bs_df)
#> Source: local data frame [6 x 4]
#> 
#>   HMP_BodySite  PHENONAME ArbitraryStatistic Entropy_bits
#>         (fctr)      (chr)              (lgl)        (dbl)
#> 1           aa  DCMCODE_1               TRUE     1.000000
#> 2           aa DCMCODE_10               TRUE     1.000000
#> 3           aa DCMCODE_11               TRUE     1.000000
#> 4           aa DCMCODE_12               TRUE     1.000000
#> 5           aa DCMCODE_13               TRUE     1.000000
#> 6           aa DCMCODE_15               TRUE     0.940286
head(ff_stats_bss_df)
#> Source: local data frame [6 x 4]
#> 
#>   HMP_BodySubsite    FUNKID Entropy_bits PercentPresent
#>            (fctr)    (fctr)        (dbl)          (dbl)
#> 1              aa  X6_44101    0.9672948      0.6486486
#> 2              aa  X5_82641    0.9866555      0.3783784
#> 3              aa X6_196156    1.0000000      0.5675676
#> 4              aa  X5_31833    0.9983637      0.5945946
#> 5              aa  X5_36281    0.9866555      0.6216216
#> 6              aa X5_105345    1.0000000      0.5675676
```

and/or save them to file, perhaps somewhere in `inspect_dir`:

``` r

save_tsv(nsamples_by_bodysite_df,
         paste0(config_l$inspect_dir, 'nsamples_by_bodysite.tsv')  # put it in `inspect_dir`
         )
```

or plot them:

``` r

plot_FFentropy_by_group(ff_stats_bss_df, 'HMP_BodySubsite')
#> `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

![](example-run_files/figure-markdown_github/unnamed-chunk-9-1.png)

Fit logistic regression models and test for coefficients significantly larger than zero
---------------------------------------------------------------------------------------

Coefficients from a logistic regression were used to identify likely associations between each subject variable and presence of each FUnkSFAM. The models account for geographic location (`SITE` variable in HMP) and were fit for each body site (or separately sub-site) and for each pair of subject variable and FUnkSFAM.

To cut down on computate time, `do_glm_tests()` calls `prefilter_by_arbitrary_statistic()` to skip model fitting if there is not enough variation in the subject variable (e.g. everyone is a vegetarian).

``` r

res_by_bodysite_df = bind_rows(
    do_glm_tests(visno_data_df, 'HMP_BodySite'),
    do_glm_tests(novisno_data_df, 'HMP_BodySite')
    )
res_by_bodysubsite_df = bind_rows(
    do_glm_tests(visno_data_df, 'HMP_BodySubsite'),
    do_glm_tests(novisno_data_df, 'HMP_BodySubsite')
    )

head(res_by_bodysite_df)
#> Source: local data frame [6 x 13]
#> 
#>   HMP_BodySite    Coeff_Name    Estimate   Std_Err     z_value   p_value
#>         (fctr)         (chr)       (dbl)     (dbl)       (dbl)     (dbl)
#> 1           aa   (Intercept) -0.19028205 0.7620272 -0.24970506 0.8028154
#> 2           aa        SITEba  0.77556670 1.0946788  0.70848792 0.4786423
#> 3           aa DCMCODE_1TRUE -0.87179358 1.0575734 -0.82433384 0.4097499
#> 4           aa   (Intercept) -0.48632014 0.7779801 -0.62510615 0.5319014
#> 5           aa        SITEba -0.09896451 1.0920914 -0.09061926 0.9277951
#> 6           aa DCMCODE_1TRUE -0.16324738 1.0379651 -0.15727637 0.8750270
#> Variables not shown: PHENONAME (chr), FUNKID (chr), N_samples (int),
#>   PH_Entropy_bits (dbl), PH_ArbitraryStatistic (lgl), FF_Entropy_bits
#>   (dbl), FF_PercentPresent (dbl)
head(res_by_bodysubsite_df)
#> Source: local data frame [6 x 13]
#> 
#>   HMP_BodySubsite    Coeff_Name   Estimate  Std_Err    z_value   p_value
#>             (chr)         (chr)      (dbl)    (dbl)      (dbl)     (dbl)
#> 1              aa   (Intercept)  0.1922806 1.417450  0.1356525 0.8920960
#> 2              aa        SITEba  0.3658118 1.399613  0.2613664 0.7938099
#> 3              aa DCMCODE_1TRUE -1.5327089 1.382456 -1.1086856 0.2675658
#> 4              aa   (Intercept) -2.3965014 1.725845 -1.3885961 0.1649556
#> 5              aa        SITEba  0.5038035 1.483917  0.3395093 0.7342261
#> 6              aa DCMCODE_1TRUE  2.0599666 1.624240  1.2682650 0.2047033
#> Variables not shown: PHENONAME (chr), FUNKID (chr), N_samples (int),
#>   PH_Entropy_bits (dbl), PH_ArbitraryStatistic (lgl), FF_Entropy_bits
#>   (dbl), FF_PercentPresent (dbl)
```

Post-filtering and p-value correction
-------------------------------------

To avoid a high penalty for multiple testing, we *a priori* decide to focus on a specific subset of tests.

### Filtering functions

<table style="width:164%;">
<colgroup>
<col width="70%" />
<col width="93%" />
</colgroup>
<thead>
<tr class="header">
<th align="left">name</th>
<th align="left">removes tests for...</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left"><code>filter_results_by_FFentropy_per_group</code></td>
<td align="left">mostly present or absent FUnkSFAMs</td>
</tr>
<tr class="even">
<td align="left"><code>filter_results_by_funksfam_annotation</code></td>
<td align="left">phylogenetically narrow FUnkSFAMs or with 10 or more annotations</td>
</tr>
<tr class="odd">
<td align="left"><code>filter_results_by_per_test_arbitrary_statistic</code></td>
<td align="left">subject variable with few differeing values</td>
</tr>
<tr class="even">
<td align="left"><code>filter_results_by_interesting_DCMCODE</code></td>
<td align="left">DCMCODEs we decided not to investigate</td>
</tr>
<tr class="odd">
<td align="left"><code>filter_unneeded_tests</code></td>
<td align="left"><code>SITE</code> and <code>(Intercept)</code> coefficients</td>
</tr>
<tr class="even">
<td align="left"><code>filter_results_all</code></td>
<td align="left">all of the above</td>
</tr>
</tbody>
</table>

#### FUnkSFAM filtering

##### `filter_results_by_funksfam_annotation`

FUnkSFAMs of span 6 and 7 were selected for analysis if they had less than 10 annotations across protein family data bases. Of these, FUnkSFAMs were included only if their presence varied appreciably among samples within body sites (or sub-sites) in order to exclude FUnkSFAMs that were present or absent from a majority of samples. See `span67_fn`.

##### `filter_results_by_FFentropy_per_group`

We used entropy (`entropy` package in CRAN) to quanitfy our uncertainty about each FUnkSFAM's presence across samples. FUnkSFAMs in the top quartile by entropy were considered for analysis.

#### Metadata filtering

##### `filter_results_by_per_test_arbitrary_statistic`

Subject variables were excluded if they did not vary appreciably among samples within body sites (or sub-sites). Subject variables were required to have at least two values with more than four observations. For example, we would choose to keep `DSUDIET` for analysis if:

> There are 5 samples with `DSUDIET = 1`, 6 with `DSUDIET = 2`, 3 with `DSUDIET = 3`

but exclude from analysis if

> There are 5 samples with `DSUDIET = 1`, 4 with `DSUDIET = 2`, 3 with `DSUDIET = 3`

### Adjusting p-values for false discovery

``` r

# apply all filters
res_bss_fil_df = filter_results_all(res_by_bodysubsite_df, ff_stats_bss_df,
                                    'HMP_BodySubsite', config_l
                                    )
#> Joining by: "HMP_BodySubsite"
#> Joining by: c("HMP_BodySubsite", "FUNKID")
#> Warning in inner_join_impl(x, y, by$x, by$y): joining factor and character
#> vector, coercing into character vector
#> Warning in inner_join_impl(x, y, by$x, by$y): joining factor and character
#> vector, coercing into character vector
#> Joining by: "NoX_FFID"

# calls p.adjust(p, method = 'fdr') and sorts the results table
res_bss_df = res_bss_fil_df %>% adjust_pvalues() %>% format_final_results()

head(res_bss_df)
#> Source: local data frame [5 x 13]
#> 
#>   HMP_BodySubsite    FUNKID PHENONAME Coeff_Name p_adj_fdr    p_value
#>             (chr)     (chr)     (chr)      (chr)     (dbl)      (dbl)
#> 1              aa X6_196156   DSUBFED    DSUBFED 0.3748134 0.07496268
#> 2              aa X6_196156   BMI_CAT  BMI_CATca 0.6049010 0.28946628
#> 3              aa X6_196156   DSUDIET    DSUDIET 0.6049010 0.43837389
#> 4              aa X6_196156   BMI_CAT  BMI_CATba 0.6049010 0.60400795
#> 5              aa X6_196156    MURICA MURICATRUE 0.6049010 0.60490098
#> Variables not shown: N_samples (int), FF_Entropy_bits (dbl),
#>   FF_PercentPresent (dbl), PH_Entropy_bits (dbl), Estimate (dbl), Std_Err
#>   (dbl), z_value (dbl)
```

``` r

# Create output directory, and save results
dir.create(config_l$output_dir, showWarnings = FALSE)
save_tsv(res_bss_df, paste0(config_l$output_dir, 'adjusted_results_by_bodysubsite.tsv'))
```
