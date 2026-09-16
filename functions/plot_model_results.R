#####################################################################
#
# Plotting function
#
# Author: F.Dhellemmes
# Last upd: 11.08.2026
#
#####################################################################
plotfunc <- function(posterior, preprandom, features_idx, x_lab, col, all_labels,
                     xlim=c(NA,NA), ylim=c(NA,NA), cex, beta_sd = NULL){
  # Build the list dynamically, scaling by beta_sd where available
  var <- setNames(
    lapply(seq_along(features_idx), function(i) {
      idx <- features_idx[i]
      draws <- posterior[[paste0("beta[", idx, "]")]]
      draws <- unlist(draws)
      if (!is.null(beta_sd) && as.character(idx) %in% names(beta_sd)) {
        draws <- draws * beta_sd[[as.character(idx)]]
      }
      draws
    }),
    paste0("var", seq_along(features_idx))
  )
  
  densities <- lapply(var, density)
  cis <- lapply(var, quantile, probs = c(0.025, 0.975))
  cis2 <- lapply(var, quantile, probs = c(0.17, 0.83))
  cex <- cex
  medians <- lapply(var, median)
  
  y_pos <- c(1:length(var))   # vertical positions
  if (is.na(ylim[1])){  ylim<-c(min(y_pos)-0.7, max(y_pos)+0.5)
  }else{ ylim<-c(ylim[1]-0.7, ylim[2]+0.5)}
  
  colors <- rep(col,length(var))
  
  if(is.na(xlim[1])){
    xlim <- range(sapply(densities, function(d) range(d$x)))
  }
  plot(NA,
       xlim = xlim,
       ylim = ylim,
       yaxt = "n",
       xlab = x_lab,
       ylab = "",
       cex.lab=cex,
       cex.axis=cex)
  
  for (i in 1:length(var)) {
    d <- densities[[i]]
    ci <- cis[[i]]
    ci2 <- cis2[[i]]
    
    y <- d$y / max(d$y) * 0.35 + y_pos[i]
    
    lines(d$x, y, col = colors[i], lwd = 1)
    lines(c(min(d$x), max(d$x)), y[1:2], col = colors[i], lwd = 1)
    lines(c(ci[1], ci[2]), y[1:2], col = colors[i], lwd = 5, lend=1)
    
    points(medians[[i]], y_pos[i], pch=19, col=colors[i], cex=1.6)
    
    polygon(
      c(d$x, rev(d$x)),
      c(y, rep(y_pos[i], length(d$x))),
      col = adjustcolor(colors[i], alpha.f = 0.2),
      border = NA
    )
    
    rand_y<-y_pos[i]-(preprandom[which(preprandom$feature == features_idx[i]),"group"]/20)
    
    points((preprandom[which(preprandom$feature == features_idx[i]),"mean_effect"]+
              preprandom[which(preprandom$feature == features_idx[i]),"mean_effect_offset"]), rand_y, pch=19, cex=1, col=scales::alpha("darkgrey",.5))
    
    for (j in 1:length(preprandom[which(preprandom$feature == features_idx[i]),"mean_effect"])){
      segments(
        x0 = (preprandom[which(preprandom$feature == features_idx[i]),"mean_effect"][j]+
                preprandom[which(preprandom$feature == features_idx[i]),"lower_effect_offset"][j]),
        y0 = rand_y[j],
        x1 = (preprandom[which(preprandom$feature == features_idx[i]),"mean_effect"][j]+
                preprandom[which(preprandom$feature == features_idx[i]),"upper_effect_offset"][j]),
        y1 = rand_y[j],
        col=scales::alpha("darkgrey",.5),
        lwd = 1
      )}
  }
  
  axis(2,
       at = y_pos,
       labels = all_labels,
       las = 1, cex.axis=cex)
  
  abline(v=0, las=2, lty=2, col="darkgrey")
}