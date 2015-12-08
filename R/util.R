# Utility functions

#' Negate shortcut
#'
#' Opposite of %in%
`%ni%` = Negate(`%in%`)

#' Cat a message to stderr
warn = function(message){
    cat(message, file = stderr())
}
