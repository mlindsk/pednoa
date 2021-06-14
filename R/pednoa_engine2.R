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
    g <- unlist(gts_pairs[k, ])
    g1 <- unique(.split_genotype(g[1]))
    g2 <- unique(.split_genotype(g[2]))
    child_genotypes <- apply(
      child <- expand.grid(g1, g2, stringsAsFactors = FALSE), 1L,
      function(x) {
        paste(as.character(sort(as.integer(x))), collapse = "/")
      }
    )
    ucg <- unique(child_genotypes)
    if (length(ucg) == 3L) {
      list(
        c(lookup[[g[1]]], lookup[[g[2]]], lookup[[ucg[1]]], 4L),
        c(lookup[[g[1]]], lookup[[g[2]]], lookup[[ucg[2]]], 2L),
        c(lookup[[g[1]]], lookup[[g[2]]], lookup[[ucg[3]]], 4L)
      )
    } else {
      lapply(ucg, function(cg) {
        c(lookup[[g[1]]], lookup[[g[2]]], lookup[[cg]], length(ucg))
      })
    }
  })

  x  <- matrix(NA, nrow = 3, ncol = nzc(na))
  v  <- vector("numeric", length = nzc(na))
  dn <- structure(list(gts, gts, gts), names = names)

  k <- 1L
  for (cell in all) {
    for (c_ in cell) {
      x[, k] <- c_[1:3]
      v[k]   <- c_[4]
      k <- k + 1L
    }
  }

  sparta::sparta_struct(x, 1/v, dn)  
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

# # na <- 5
# # nh <- choose(na, 2)
# # nzc(na)
# # gt <- make_genotypes(na)


# # all <- apply(gt_pairs, 1L, function(g) {
# #   g1 <- unique(.split_genotype(g[1]))
# #   g2 <- unique(.split_genotype(g[2]))
# #   child_genotypes <- apply(
# #     child <- expand.grid(g1, g2, stringsAsFactors = FALSE), 1L,
# #     function(x) {
# #       paste(as.character(sort(as.integer(x))), collapse = "/")
# #     }
# #   )
# #   ucg <- unique(child_genotypes)
# #   if (length(ucg) == 3L) {
# #     list(
# #       c(lookup[[g[1]]], lookup[[g[2]]], lookup[[ucg[1]]], 4L),
# #       c(lookup[[g[1]]], lookup[[g[2]]], lookup[[ucg[2]]], 2L),
# #       c(lookup[[g[1]]], lookup[[g[2]]], lookup[[ucg[3]]], 4L)
# #     )
# #   } else {
# #     lapply(ucg, function(cg) {
# #       c(lookup[[g[1]]], lookup[[g[2]]], lookup[[cg]], length(ucg))
# #     })    
# #   }
# # })


# # l <- all |> .map_int(length)
# # sum(l)

# # which(l == 1) |> length() == na^2
# # which(l == 2) |> length() == 2 * na * nh
# # which(l == 3) |> length() == nh
# # which(l == 4) |> length() == nh * (nh - 1)
