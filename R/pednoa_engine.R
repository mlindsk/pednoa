# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                    HELPERS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.split_genotype <- function(x) strsplit(x, "/")[[1]]

.unique_genotypes_from_parents <- function(p1_genotype, p2_genotype) {
  combs. <- expand.grid(
    .split_genotype(p1_genotype),
    .split_genotype(p2_genotype),
    stringsAsFactors = FALSE
  )
  unique(apply(combs., 1L, function(x) paste0(sort(as.numeric(x)), collapse = "/")))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#               CPT CONSTRUCTORS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## af <- c(.1,.2,.7)
## ua = .map_chr(1:length(af), as.character)  

cpt_founder <- function(name, af) {
  af_names <- seq_along(af)
  het_comb <- utils::combn(af_names, 2)
  hw_prob  <- c(af^2, apply(het_comb, 2L, function(x) af[x[1]] * af[x[2]]) * 2)
  het_names <- apply(het_comb, 2L, paste, collapse = "/")
  hom_names <- paste(af_names, af_names, sep = "/")
  genotypes <- c(hom_names, het_names)
  array(
    data = hw_prob,
    dim = length(genotypes),
    dimnames = structure(list(genotypes), names = name)
  )
}

## cpt_founder2(c("Founder"), af)

cpt_child_parents <- function(names, gts){
  # The child name must be the first element in names
  # gts: as the names returned from cpt_fouders

  ngts <- length(gts)
  A <- array(0L,
    dim = rep(ngts, 3L),
    dimnames = structure(list(gts, gts, gts), names = names)
  )

  dim_A12 <- dimnames(A)[1:2]
  
  for(g2 in gts) {
    A_g2 <- matrix(0L, ngts, ngts)
    dimnames(A_g2) <- dim_A12
    for (g1 in gts) {
      ugts <- .unique_genotypes_from_parents(g1, g2)
      A_g2[ugts, g1] <- 1 / length(ugts)
    }
    A[, , g2] <- A_g2
  }

  A
}

## cpt_child_parents2(c("C", "P1", "P2"), names(cpt_founder(c("Founder"), af)))

make_founder_cpts <- function(founder_vec, afreq) {
  generic_cpt_founder <- cpt_founder(0, afreq)
  lapply(founder_vec, function(x) {
    cpt <- generic_cpt_founder
    names(dimnames(cpt)) <- x
    cpt
  })
}

make_child_parents_cpts <- function(genotypes, nonfounder_list) {
  generic_cpt_child_parents <- cpt_child_parents(1:3, genotypes)
  lapply(nonfounder_list, function(x) {
    cpt <- generic_cpt_child_parents
    names(dimnames(cpt)) <- x
    cpt
  })
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#        EXTRACT NUMBER OF DIFFERENT ALLELES
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.n_unique_col <- function(x, j) {
  cell <- sparta::get_cell_name(x, j)
  length(unique(unlist(strsplit(cell, "/"))))
}

.pmf_unique_alleles <- function(j, mvec_chr, afreq) {

  mvec_clique_idx <- which(.map_lgl(j$cliques, function(x) {
    all(mvec_chr %in% x)
  }))[1]
  
  mvec_clique <- j$charge$C[[mvec_clique_idx]]
  mvec_pmf <- sparta::marg(mvec_clique, setdiff(names(mvec_clique), mvec_chr))

  # TODO: Can this be speeded up? Possibly with some C code yes (so do it?!)
  number_of_different_alleles <- apply(
    mvec_pmf,
    2L,
    function(j) .n_unique_col(mvec_pmf, j)
  )

  unique_ndifferent_alleles <- unique(number_of_different_alleles)
  sorted_unique_ndifferent_alleles <- sort(unique_ndifferent_alleles)
  
  p <- .map_dbl(sorted_unique_ndifferent_alleles, function(u) {
    u_idx <- which(number_of_different_alleles == u)
    sum(attr(mvec_pmf, "vals")[u_idx])
  })

  possible_ndiff_alleles  <- 1:length(afreq)
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
                   trace = TRUE
                   ) {


  mvec_chr <- as.character(mvec)

  founder_cpts    <- make_founder_cpts(founder_vec, afreq)
  nonfounder_cpts <- make_child_parents_cpts(names(founder_cpts[[1L]]), nonfounder_list)

  cpts <- structure(
    c(founder_cpts, nonfounder_cpts),
    names = c(founder_vec, .map_dbl(nonfounder_list, function(x) x[1]))
  )
  
  cl <- jti::cpt_list(cpts)
  rm(cpts)

  cp <- jti::compile(cl, joint_vars = mvec_chr, save_graph = TRUE, opt = "min_fill")

  e <- unlist(lapply(evidence, function(e) {
    structure(paste0(sort(e[2:3]), collapse = "/"), names = e[1])
  }))


  # NOTE: Test if evidence is correct from a pedigree point of view!
  # Just test for zero sparta? 
  j    <- jti::jt(cp, e)
  pout <- .pmf_unique_alleles(j, mvec_chr, afreq)
  
  structure(
    as.array(pout),
    dimnames = structure(
      list(names(pout)),
      names = "P(Nm)"
    )
  )
}
