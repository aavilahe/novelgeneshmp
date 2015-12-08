# Utility functions

#' Negate shortcut
#'
#' Opposite of %in%
#' @export
`%ni%` = Negate(`%in%`)

#' Cat a message to stderr
#' @export
warn = function(message){
    cat(message, file = stderr())
}
