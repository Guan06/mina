################################################################################
# Define all accessors to get and set slots

#' Setter and getter for the slot `tab`.
#' @param x The `mina` object.
#' @param value The value to set for the slot of the `mina` object `x`.
#' @examples
#' tab(maize) <- maize_asv2
#' tab(maize)[1:5, 1:5]
#' @rdname tab_accessor
#' @export

setGeneric("tab<-", function(x, value) standardGeneric("tab<-"))

#' @rdname tab_accessor
#' @exportMethod tab<-

setMethod("tab<-", "mina", function(x, value) {
              x@tab <- value
              x
})

#' @rdname tab_accessor
#' @export

setGeneric("tab", function(x) standardGeneric("tab"))

#' @rdname tab_accessor
#' @exportMethod tab

setMethod("tab", "mina", function(x) x@tab)

#' Setter and getter for the slot `des`, which is the description and meta data
#' of rows in `tab`.
#' @param x The `mina` object.
#' @param value The value to set for the slot of the `mina` object `x`.
#' @examples
#' des(maize) <- maize_des2
#' head(des(maize))
#' @rdname des_accessor
#' @export

setGeneric("des<-", function(x, value) standardGeneric("des<-"))

#' @rdname des_accessor
#' @exportMethod des<-

setMethod("des<-", "mina", function(x, value) {
              x@des <- value
              x
})

#' @rdname des_accessor
#' @export
setGeneric("des", function(x) standardGeneric("des"))

#' @rdname des_accessor
#' @exportMethod des

setMethod("des", "mina", function(x) x@des)

#' Setter and getters for the slot `norm`, normalized `tab` matrix.
#' @param x The `mina` object.
#' @param value The value to set for the slot of the `mina` object `x`.
#' @examples
#' norm(maize) <- norm_tab(maize_asv2, method = "total")
#' norm(maize)[1:5, 1:5]
#' @rdname norm_accessor
#' @export

setGeneric("norm<-", function(x, value) standardGeneric("norm<-"))

#' @rdname norm_accessor
#' @exportMethod norm<-

setMethod("norm<-", "mina", function(x, value) {
              x@norm <- value
              x
})

#' @rdname norm_accessor
#' @export

setGeneric("norm", function(x) standardGeneric("norm"))

#' @rdname norm_accessor
#' @exportMethod norm

setMethod("norm", "mina", function(x) x@norm)

#' Setter for the slot `adj` and `adj_sig`, the adjacency matrix of
#' `norm` and corresponding significant value matrix with `sig` is `TRUE`.
#' @param x The `mina` object.
#' @param value The value to set for the slot of the `mina` object `x`.
#' @examples
#' norm(maize) <- norm_tab(maize_asv2, method = "total")
#' maize <- adj(maize, method = "spearman", sig = TRUE)
#' adj_sig(maize)[1:5, 1:5]
#' @rdname adj_accessor
#' @keywords internal

setGeneric(".adj<-", function(x, value) standardGeneric(".adj<-"))

#' @rdname adj_accessor
#' @keywords internal

setMethod(".adj<-", "mina", function(x, value) {
              x@adj <- value
              x
})

#' Get the slot `adj`.
#' @param x The `mina` object.
#' @rdname adj_accessor
#' @keywords internal

.adj <- function(x) x@adj

#' @rdname adj_accessor

setGeneric("adj_sig<-", function(x, value) standardGeneric("adj_sig<-"))

#' @rdname adj_accessor
#' @keywords internal

setMethod("adj_sig<-", "mina", function(x, value) {
              x@adj_sig <- value
              x
})

#' Get the slot `adj_sig`.
#' @param x The `mina` object.
#' @keywords internal

adj_sig <- function(x) x@adj_sig

#' Setter for the slot `cls`.
#' @param x The `mina` object.
#' @param value The value to set for the slot of the `mina` object `x`.
#' @rdname cls_accessor
#' @keywords internal

setGeneric("cls<-", function(x, value) standardGeneric("cls<-"))

#' @rdname cls_accessor
#' @keywords inter

setMethod("cls<-", "mina", function(x, value) {
              x@cls <- value
              x
})

#' Get the slot `cls`.
#' @param x The `mina` object.
#' @keywords internal

cls <- function(x) x@cls

#' Setter for the slot `cls_tab`.
#' @param x The `mina` object.
#' @param value The value to set for the slot of the `mina` object `x`.
#' @rdname cls_tab_accessor
#' @keywords internal

setGeneric("cls_tab<-", function(x, value) standardGeneric("cls_tab<-"))

#' @rdname cls_tab_accessor
#' @keywords internal

setMethod("cls_tab<-", "mina", function(x, value) {
              x@cls_tab <- value
              x
})

#' Get the slot `cls_tab`.
#' @param x The `mina` object.
#' @keywords internal

cls_tab <- function(x) x@cls_tab

#' Setter and getter for the slot `dis`.
#' @param x The `mina` object.
#' @param value The value to set for the slot of the `mina` object `x`.
#' @rdname dis_accessor
#' @export

setGeneric("dis<-", function(x, value) standardGeneric("dis<-"))

#' @rdname dis_accessor
#' @export

setGeneric("dis", function(x) standardGeneric("dis"))

#' @rdname dis_accessor
#' @exportMethod dis<-

setMethod("dis<-", "mina", function(x, value) {
              x@dis <- value
              x
})

#' Get the slot `dis`
#' @param x The `mina` object.
#' @rdname dis_accessor
#' @exportMethod dis

setMethod("dis", "mina", function(x) x@dis)

#' Setter and getter for the slot `dmr`.
#' @param x The `mina` object.
#' @param value The value to set for the slot of the `mina` object `x`.
#' @rdname dmr_accessor

setGeneric(".dmr<-", function(x, value) standardGeneric(".dmr<-"))

#' @rdname dmr_accessor

setGeneric(".dmr", function(x) standardGeneric(".dmr"))

#' @rdname dmr_accessor
#' @keywords internal

setMethod(".dmr<-", "mina", function(x, value) {
              x@dmr <- value
              x
})

#' @rdname dmr_accessor
#' @keywords internal

setMethod(".dmr", "mina", function(x) x@dmr)

#' Setter and getter for the slot `multi` and `perm`.
#' @param x The `mina` object.
#' @param value The value to set for the slot of the `mina` object `x`.
#' @rdname bs_pm_accessor

setGeneric("multi<-", function(x, value) standardGeneric("multi<-"))

#' @rdname bs_pm_accessor

setGeneric("multi", function(x) standardGeneric("multi"))

#' @rdname bs_pm_accessor
#' @keywords internal

setMethod("multi<-", "mina", function(x, value) {
              x@multi <- value
              x
})

#' @rdname bs_pm_accessor
#' @keywords internal

setMethod("multi", "mina", function(x) x@multi)

#' @rdname bs_pm_accessor

setGeneric("perm<-", function(x, value) standardGeneric("perm<-"))

#' @rdname bs_pm_accessor

setGeneric("perm", function(x) standardGeneric("perm"))

#' @rdname bs_pm_accessor
#' @keywords internal

setMethod("perm<-", "mina", function(x, value) {
              x@perm <- value
              x
})

#' @rdname bs_pm_accessor
#' @keywords internal

setMethod("perm", "mina", function(x) x@perm)

#' Setter for the slots 'dis_bs', `dis_pm` and `dis_stat`.
#' @param x The `mina` object.
#' @param value The value to set for the slot of the `mina` object `x`.
#' @rdname net_dis_accessor
#' @keywords internal

setGeneric("dis_bs<-", function(x, value) standardGeneric("dis_bs<-"))

#' @rdname net_dis_accessor

setGeneric("dis_pm<-", function(x, value) standardGeneric("dis_pm<-"))

#' @rdname net_dis_accessor

setGeneric("dis_stat<-", function(x, value) standardGeneric("dis_stat<-"))

#' @rdname net_dis_accessor
#' @keywords internal

setMethod("dis_bs<-", "mina", function(x, value) {
              x@dis_bs <- value
              x
})

#' @rdname net_dis_accessor
#' @keywords internal

setMethod("dis_pm<-", "mina", function(x, value) {
              x@dis_pm <- value
              x
})

#' @rdname net_dis_accessor
#' @keywords internal

setMethod("dis_stat<-", "mina", function(x, value) {
              x@dis_stat <- value
              x
})

#' Getter for the slots `dis_bs`, `dis_pm` and `dis_stat`.
#' @param x The `mina` object.
#' @export

dis_bs <- function(x) x@dis_bs
dis_pm <- function(x) x@dis_pm
dis_stat <- function(x) x@dis_stat
