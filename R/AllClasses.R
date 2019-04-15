#' Class "mina" includes the quantitative table and descriptive table.
#'
#'
#' @examples
#' new(mina)
#' @export
setClass('mina',
         representation(Tab = "matrix",
                   desTab = "data.frame",
                   normTab = "matrix",
                   dmr = "list"
         )
         )
