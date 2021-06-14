# make_genotypes <- function(na) {
#   het_comb  <- utils::combn(1:na, 2)
#   het_names <- apply(het_comb, 2L, paste, collapse = "/")
#   hom_names <- .map_chr(1:na, function(a) paste(a, a, sep = "/"))
#   c(hom_names, het_names)
# }

# nzc <- function(na) {
#   nh <- choose(na, 2)
#   ng <- nh + na
#   na^2 +  4*na*nh + 3*nh + 4*nh*(nh - 1)
# }

# na <- 5
# nh <- choose(na, 2)
# nzc(na)
# gt <- make_genotypes(na)

# # make a lookup table that maps genotypes to an index
# lookup <- new.env()
# for (i in seq_along(gt)) lookup[[gt[i]]] <- i
# gt_pairs <- expand.grid(gt, gt, stringsAsFactors = FALSE)

# # future::plan(multisession, workers = 3)
# # future::plan(sequential)
# # TODO: Just mclapply over 1:nrow(gt_pairs)

# all <- apply(gt_pairs, 1L, function(g) {
#   g1 <- unique(.split_genotype(g[1]))
#   g2 <- unique(.split_genotype(g[2]))
#   child_genotypes <- apply(
#     child <- expand.grid(g1, g2, stringsAsFactors = FALSE), 1L,
#     function(x) {
#       paste(as.character(sort(as.integer(x))), collapse = "/")
#     }
#   )
#   ucg <- unique(child_genotypes)
#   # TODO: Repair for length(ucg) == 3 (two with 1/4 and one with 1/2)
#   lapply(ucg, function(cg) {
#     c(lookup[[g[1]]], lookup[[g[2]]], lookup[[cg]], length(ucg))
#   })
# })


# l <- all |> .map_int(length)
# sum(l)

# # which(l == 1) |> length() == na^2
# # which(l == 2) |> length() == 2 * na * nh
# # which(l == 3) |> length() == nh
# # which(l == 4) |> length() == nh * (nh - 1)

# x  <- matrix(NA, nrow = 3, ncol = nzc(na))
# v  <- vector("numeric", length = nzc(na))
# dn <- list(P1 = gt, P2 = gt, C = gt)

# # print(TRUE)

# k <- 1L
# for (cell in all) {
#   for (c_ in cell) {
#     x[, k] <- c_[1:3]
#     v[k]   <- c_[4]
#     k <- k + 1L
#   }
# }

# any(is.na(v))

# sparse_pednoa_cpt <- sparta::sparta_struct(x, 1/v, dn)
