## MAPS
.map_chr <- function(x, fun, ...) vapply(X = x, FUN = fun, FUN.VALUE = character(1), ...)
.map_int <- function(x, fun, ...) vapply(X = x, FUN = fun, FUN.VALUE = integer(1), ...)
.map_dbl <- function(x, fun, ...) vapply(X = x, FUN = fun, FUN.VALUE = numeric(1), ...)
.map_lgl <- function(x, fun, ...) vapply(X = x, FUN = fun, FUN.VALUE = logical(1), ...)
neq_null <- function(x) !is.null(x)
'%ni%'   <- Negate('%in%')
eq_empt_lst <- function(x) inherits(x, "list") && length(x) == 0L

## ASSERTERS
sum_to_one_ <- function(afreq_list) {
  sum(.map_dbl(afreq_list, function(x) sum(x) > 0.999999)) == length(afreq_list)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                    HELPERS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.split_genotype <- function(x) strsplit(x, "/")[[1]]

is_whole_number <- function(x) identical(round(x), x)

is_child_parent <- function(x, fv) {
  length(x) == 3L &&
    is_whole_number(x) &&
    length(unique(x)) == length(x) &&
    x[1] %ni% fv
}

extract_evidence <- function(evidence_list, nloci) {
  ids <- lapply(1:nloci, function(l) {
    lapply(evidence_list, function(e) {
      id <- e[[1]]
      loci <- e[[2]]
      if (l %ni% loci) return(NULL)
      loci_l_idx <- match(l, loci)
      alleles_l <- e[[3]][[loci_l_idx]]
      c(id, alleles_l)
    })
  })
  lapply(ids, function(x) {
    f <- Filter(function(z) !is.null(z), x)
    if (eq_empt_lst(f)) NULL else f
  })
}
