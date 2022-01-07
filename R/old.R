# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#               CPT CONSTRUCTORS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# af <- c(.1,.2,.7, .1, .1)
# ua = .map_chr(1:length(af), as.character)  

# cpt_founder <- function(name, af) {
#   af_names <- seq_along(af)
#   het_comb <- utils::combn(af_names, 2)
#   hw_prob  <- c(af^2, apply(het_comb, 2L, function(x) af[x[1]] * af[x[2]]) * 2)
#   het_names <- apply(het_comb, 2L, paste, collapse = "/")
#   hom_names <- paste(af_names, af_names, sep = "/")
#   genotypes <- c(hom_names, het_names)
#   array(
#     data = hw_prob,
#     dim = length(genotypes),
#     dimnames = structure(list(genotypes), names = name)
#   )
# }

# # cpt_founder(c("Founder"), c(.1, .2, .3, .4))

# cpt_child_parents <- function(names, gts){
#   # The child name must be the first element in names
#   # gts: as the names returned from cpt_fouders

#   ngts <- length(gts)
#   A <- array(0L,
#     dim = rep(ngts, 3L),
#     dimnames = structure(list(gts, gts, gts), names = names)
#   )

#   dim_A12 <- dimnames(A)[1:2]
  
#   for(g2 in gts) {
#     A_g2 <- matrix(0L, ngts, ngts)
#     dimnames(A_g2) <- dim_A12
#     for (g1 in gts) {
#       ugts <- .unique_genotypes_from_parents(g1, g2)
#       A_g2[ugts, g1] <- 1 / length(ugts)
#     }
#     A[, , g2] <- A_g2
#   }

#   A
# }

# make_founder_cpts <- function(founder_vec, afreq) {
#   generic_cpt_founder <- cpt_founder(0, afreq)
#   lapply(founder_vec, function(x) {
#     cpt <- generic_cpt_founder
#     names(dimnames(cpt)) <- x
#     cpt
#   })
# }

# make_child_parents_cpts <- function(genotypes, nonfounder_list) {
#   generic_cpt_child_parents <- cpt_child_parents(1:3, genotypes)
#   lapply(nonfounder_list, function(x) {
#     cpt <- generic_cpt_child_parents
#     names(dimnames(cpt)) <- x
#     cpt
#   })
# }

# .unique_genotypes_from_parents <- function(p1_genotype, p2_genotype) {
#   combs. <- expand.grid(
#     .split_genotype(p1_genotype),
#     .split_genotype(p2_genotype),
#     stringsAsFactors = FALSE
#   )
#   unique(apply(combs., 1L, function(x) paste0(sort(as.numeric(x)), collapse = "/")))
# }

