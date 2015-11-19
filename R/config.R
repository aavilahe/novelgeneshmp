
# config documentation
load_config = function(json_fn){
    # loads a config
    config_l = jsonlite::fromJSON(file(json_fn))
    return(config_l)
}


