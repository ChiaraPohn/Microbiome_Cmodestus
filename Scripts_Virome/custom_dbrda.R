#' Functions are adapted or copied (because they were not exported) from https://github.com/raeslab/RLdbRDA/tree/main
custom_rldbrda <- function(distmat, meta, p_cutoff = 0.05) {
  r2 <- get_r2(distmat, meta)

  sign_r2 <- rownames(r2[which(r2$padj < p_cutoff), ]) # selects only variables significant in step 1

  if (length(sign_r2) < 1) {
    message("No significant features found!")
    return(NULL)
  }

  cumul <- get_cumul(distmat, meta[, sign_r2, drop = F])

  out <- combine_data(r2, cumul)

  return(out)
}

get_r2_single <- function(distmat, meta, feature) {
  capsc <- vegan::capscale(distmat ~ meta[, feature], na.action = na.omit)
  an <- vegan::anova.cca(capsc, permutations = 9999)

  Fa <- an["F"][[1]][[1]]
  r2 <- vegan::RsquareAdj(capsc)[[1]]
  r2adj <- vegan::RsquareAdj(capsc)[[2]]
  N <- nrow(na.exclude(meta[, feature, drop = FALSE]))
  pval <- an["Pr(>F)"][[1]][[1]]

  output <- cbind(feature, Fa, r2, r2adj, N, pval)

  return(output)
}

get_r2 <- function(distmat, meta) {
  features <- colnames(meta)

  all <- c()

  for (feature in features) {
    es <- get_r2_single(distmat, meta, feature)
    all <- rbind(all, es)
  }

  all <- data.frame(all)
  all$padj <- p.adjust(all$pval, method = "BH")

  rownames(all) <- all$feature

  return(all)
}

get_cumul <- function(distmat, meta) {
  mod0 <- vegan::capscale(distmat ~ 1) # H0: unconstrained ordination
  mod1 <- vegan::capscale(distmat ~ ., data = meta) # H1: full constrained ordination, all metadata


  attach(meta)

  step.res <- vegan::ordiR2step(mod0, scope = formula(mod1), data = meta, direction = "forward", Pin = 1, R2scope = TRUE, pstep = 100, permutations = 9999, trace = F) # forward stepwise dbRDA
  res <- step.res$anova


  row.names(res) <- gsub(pattern = "\\+ ", "", row.names(res))
  colnames(res) <- gsub(pattern = "Pr\\(>F\\)", "pval", colnames(res)) # replace column name
  colnames(res) <- paste0("RDAcumul_", colnames(res))
  res[, "RDAcumul_N"] <- nrow(meta)

  detach(meta)


  return(res)
}

combine_data <- function(r2, cumul) {
  all <- data.frame(merge(r2, cumul, by = "row.names", all = T), row.names = 1)
  all <- all[order(all$r2, decreasing = TRUE), ]
  all <- all[order(all$RDAcumul_R2.adj), ]

  return(all)
}
