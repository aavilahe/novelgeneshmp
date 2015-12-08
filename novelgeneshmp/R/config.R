# Config file loading functions

#' Loads config file
#'
#' NOTE: requires \code{jsonlite} package
#'
#' @param json_fn
#' @return A list of input/output file paths and parameters
#' @export
load_config = function(json_fn){
    # loads a config
    config_l = jsonlite::fromJSON(file(json_fn))
    return(config_l)
}


