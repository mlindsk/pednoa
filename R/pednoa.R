#' Pedigree Number of different Alleles
#'
#' Find the probability distribution of m different pedigree members having
#' different number of alleles
#'
#' @param afreq_list A list with allele frequencies. One for each DNA loci
#' @param mvec Integer vector with the index of the pedigree members of interest
#' @param founder_vec Integer vector with the index of all founder
#' @param nonfounder_list List of three-dimensional index vectors of all
#' nonfounders on the form \code{(child, paternal_index, maternal_index)}
#' @param evidence_list A list of evidence. See the examples for details.
#' @param locus_wise Logical. \code{FALSE} if the pmfs should be convoluted
#' @param ncores Integer. The number of cores to be used in connection with
#' parallelizing some of the internal code.
#' @details There is no check whether the evidence is in agreement with
#' the pedigree. The user is responsible for entering correct evidence.
#' @examples
#'
#' # Allele frequencies:
#' # -------------------
#' 
#' p1 <- rbeta(6, 3, 6)
#' p1 <- p1 / sum(p1)
#' 
#' p2 <- rbeta(6, 3, 6)
#' p2 <- p2 / sum(p2)
#' 
#' # Evidence:
#' # -------------------
#'
#' e <- list(
#'   list(3, c(1, 2), list(c(2, 2), c(1, 3))),
#'   list(1, 2, list(c(2, 3)))
#' )
#'
#' # Individual 3 has genotype (2,2) on loci 1 and (1,3) on loci two.
#' # Individual 1 has genotype (2,2) on loci 2
#'
#'
#' # Consider the pedigree:
#' # ----------------------
#' 
#' #  1-------2
#' #      |
#' #  ---------
#' #  |       |
#' #  3       4-------6
#' #              |
#' #              5
#'
#' # We let 1 and 5 be the individuals of interest. That is, we calculate the
#' # probability of seeing different alleles over two locis for 1 and 5 given
#' # the evidence.
#' 
#' p <- pednoa(
#'   list(p1, p2),
#'   c(1, 5),
#'   c(1, 2, 6),
#'   list(c(3, 1, 2), c(4, 1, 2), c(5, 4, 6)),
#'   e,
#'   locus_wise = FALSE
#' )
#'
#' p
#' 
#' @export
pednoa <- function(
                   afreq_list,
                   mvec,      
                   founder_vec,
                   nonfounder_list,
                   evidence_list = NULL,
                   locus_wise = TRUE,
                   ncores = 1
                   ) {

  if (!is.list(afreq_list)) stop("afreq_list must be a list even if it contains a single element")
  
  if (!sum_to_one_(afreq_list)) stop("Some alle frequencies do not sum to one", call. = FALSE)
  
  if (!is_whole_number(mvec)) {
    stop("mvec must be a vector of integers")
  }

  if (!is_whole_number(founder_vec)) {
    stop("founder_vec must be a vector of integers")
  }

  if (!all(.map_lgl(nonfounder_list, is_child_parent, founder_vec))) {
    stop("some elements in nonfounder_list are not on correct form")
  }
  
  nloci <- length(afreq_list)
  e <- extract_evidence(evidence_list, nloci)

  pmfs <- lapply(seq_along(afreq_list), function(l) {
    afreq_l <- afreq_list[[l]]
    el <- e[[l]]
    pednoa_engine(afreq_l, mvec, founder_vec, nonfounder_list, el, ncores)
  })

  if (nloci == 1L) return(pmfs[[1]])
  if (locus_wise) return(pmfs)
  
  nmvec <- length(mvec)
  attr(pmfs[[1]], "ell") <- 2L

  convolute <- function(p1, p2) {
    ell <- attr(p1, "ell")
    n_current_max <- min(as.integer(names(p1[length(p1)])) + length(p2), ell* 2*nmvec)
    ns <- ell:n_current_max
    p3 <- structure(vector("double", length = length(ns)), names = ns)

    pairs <- expand.grid(as.integer(names(p1)), as.integer(names(p2)))
    pairs$sum_ <- pairs[, 1] + pairs[, 2]
    pairs[, 1] <- as.character(pairs[, 1])
    pairs[, 2] <- as.character(pairs[, 2])

    N_counter <- 1L
    for (N in ns) {
      idx_N <- which(pairs[, "sum_"] == N)
      pairs_N <- pairs[idx_N, , drop = FALSE]
      p3[N_counter] <- sum(apply(pairs_N, 1L, function(x) {
        p1[x[1]] * p2[x[2]]
      }))
      N_counter <- N_counter + 1L
    }
    
    attr(p3, "ell") <- ell + 1L
    p3
  }

  out <- as.array(Reduce(convolute, pmfs))
  attr(out, "ell") <- NULL
  structure(out, dimnames = structure(list(names(out)), names = "P(Nm)"))
}

# p <- pednoa(
#   list(p1, p2, p3, p4),
#   c(1, 5),
#   c(1, 2, 6),
#   list(c(3, 1, 2), c(4, 1, 2), c(5, 4, 6))
# )
