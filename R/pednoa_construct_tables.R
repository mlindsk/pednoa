nzc <- function(na) {
  nh <- choose(na, 2)
  ng <- nh + na
  na^2 +  4*na*nh + 3*nh + 4*nh*(nh - 1)
}

make_genotypes <- function(na) {
  het_comb  <- utils::combn(1:na, 2)
  het_names <- apply(het_comb, 2L, paste, collapse = "/")
  hom_names <- .map_chr(1:na, function(a) paste(a, a, sep = "/"))
  c(hom_names, het_names)
}

cpt_founder_sparse <- function(name, af) {  
  af_names <- seq_along(af)
  het_comb <- utils::combn(af_names, 2)
  hw_prob  <- c(af^2, apply(het_comb, 2L, function(x) af[x[1]] * af[x[2]]) * 2)
  gts      <- make_genotypes(length(af))
  x        <- matrix(1:length(gts), nrow = 1)
  dn       <- structure(list(gts), names = name)
  sparta::sparta_struct(x, hw_prob, dn)
}

cpt_child_parents_sparse <- function(names, na, ncores = 1) {
  
  gts <- make_genotypes(na)
  
  # make a lookup table that maps genotypes to an index
  lookup <- new.env()
  for (i in seq_along(gts)) lookup[[gts[i]]] <- i
  gts_pairs <- expand.grid(gts, gts, stringsAsFactors = FALSE)

  all <- parallel::mclapply(mc.cores = ncores, X = 1:nrow(gts_pairs), FUN = function(k) {

    g  <- unlist(gts_pairs[k, ])
    g1 <- .split_genotype(g[1])
    g2 <- .split_genotype(g[2])
    
    child_genotypes <- c(
      paste(sort(as.integer(c(g1[1], g2[1]))), collapse = "/"),
      paste(sort(as.integer(c(g1[1], g2[2]))), collapse = "/"),
      paste(sort(as.integer(c(g1[2], g2[1]))), collapse = "/"),
      paste(sort(as.integer(c(g1[2], g2[2]))), collapse = "/"))
    
    child_prob <- table(child_genotypes) / sum(table(child_genotypes))
    child_names_int <- as.integer(sapply(names(child_prob), function(x) lookup[[x]]))

    par_idx <- c(lookup[[g[1]]], lookup[[g[2]]])

    lapply(seq_along(child_names_int), function(i) {
      list(index = c(child_names_int[i], par_idx), val = unname(child_prob[i]))
    })
    
  })
  
  x  <- matrix(NA, nrow = 3, ncol = nzc(na))
  v  <- vector("numeric", length = nzc(na))
  dn <- structure(list(gts, gts, gts), names = names)

  k <- 1
  for (a in all) for (e in a) {
    x[, k] <- e$index
    v[k] <- e$val
    k <- k + 1L
  }
  
  sparta::sparta_struct(x, v, dn)  
}

make_founder_cpts_sparse <- function(founder_vec, afreq) {
  generic_cpt_founder <- cpt_founder_sparse(0, afreq)
  lapply(founder_vec, function(x) {
    cpt <- generic_cpt_founder
    names(attr(cpt, "dim_names")) <- x
    cpt
  })
}

make_child_parents_cpts_sparse <- function(nonfounder_list, na, ncores = 1) {
  generic_cpt_child_parents <- cpt_child_parents_sparse(1:3, na, ncores)
  lapply(nonfounder_list, function(x) {
    cpt <- generic_cpt_child_parents
    names(attr(cpt, "dim_names")) <- x
    cpt
  })
}


