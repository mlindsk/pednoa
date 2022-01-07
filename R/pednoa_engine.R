n_unique_col <- function(x, j) {
  cell <- sparta::get_cell_name(x, j)
  as.integer(length(unique(unlist(strsplit(cell, "/")))))
}

pmf_unique_alleles <- function(j, mvec_chr, afreq, ncores = 1) {

  mvec_clique_idx <- as.integer(substr(attr(j, "clique_root"), 2, 2))
  mvec_clique     <- j$charge$C[[mvec_clique_idx]]
  mvec_pmf        <- sparta::marg(mvec_clique, setdiff(names(mvec_clique), mvec_chr))

  number_of_different_alleles <- unlist(parallel::mclapply(
    mc.cores = ncores, X = 1:ncol(mvec_pmf), FUN = function(k) {
    cell <- mvec_pmf[, k, drop = TRUE]
    n_unique_col(mvec_pmf, cell)
  }))
  
  
  unique_ndifferent_alleles <- unique(number_of_different_alleles)
  sorted_unique_ndifferent_alleles <- sort(unique_ndifferent_alleles)
  
  p <- .map_dbl(sorted_unique_ndifferent_alleles, function(u) {
    u_idx <- which(number_of_different_alleles == u)
    sum(attr(mvec_pmf, "vals")[u_idx])
  })

  possible_ndiff_alleles <- 1:length(afreq)
  zeroe_names <- which(possible_ndiff_alleles %ni% sorted_unique_ndifferent_alleles)
  pout <- structure(rep(0, length(afreq)), names = 1:length(afreq))
  pout[sorted_unique_ndifferent_alleles] <- p
  pout
}


pednoa_engine <- function(
                   afreq,
                   mvec,             # vector of rascals indices
                   founder_vec,      # vector of founder indices
                   nonfounder_list,  # list of 3-dim vectors with indices (child, parents)
                   evidence = NULL,
                   ncores   = 1L
                   ) {

  mvec_chr        <- as.character(mvec)
  founder_cpts <- make_founder_cpts_sparse(founder_vec, afreq)
  nonfounder_cpts <- make_child_parents_cpts_sparse(nonfounder_list, length(afreq), ncores)

  cpts <- structure(
    c(founder_cpts, nonfounder_cpts),
    names = c(founder_vec, .map_chr(nonfounder_cpts, function(x) names(x)[1]))
  )

  cl <- jti::cpt_list(cpts)
  rm(cpts)

  e <- unlist(lapply(evidence, function(e) {
    structure(paste0(sort(e[2:3]), collapse = "/"), names = e[1])
  }))

  cp   <- jti::compile(cl, evidence = e, joint_vars = mvec_chr)
  j    <- jti::jt(cp, evidence = e, propagate = "collect")
  pout <- pmf_unique_alleles(j, mvec_chr, afreq, ncores)
  
  structure(
    as.array(pout),
    dimnames = structure(
      list(names(pout)),
      names = "P(Nm)"
    )
  )
}

