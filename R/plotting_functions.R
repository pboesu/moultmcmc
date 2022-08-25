### pairwise plot of marginal densities
### code following a stackexchange example

#' Pairwise plots of posterior draws
#'
#' Plots pairwise scatter/density plots of posterior draws
#'
#' @param x a moultmcmc object
#' @param y an optional second moultmcmc object
#' @param pars 	An optional character vector of parameter names. If pars is not specified then the default is to use the first 4 parameters of x.
#' @param names optional character vector of model names
#' @param scatter logical, draw scatterplot of posterior samples if overlay = FALSE, else draw contours
#' @param overlay logical, plot scatterplots and density contours of both models in each triangle, else draw one model per triangle
#' @param ... further arguments to plot.default (the call that draws the scatter/contour plot)
#' @importFrom MASS kde2d
#' @import RColorBrewer
#' @importFrom graphics abline contour hist layout lines par plot plot.default points text title strwidth legend
#' @importFrom stats density
#' @importFrom rstan extract
#' @export
compairs_plot <- function(x, y = NULL, pars = NULL, scatter = TRUE, overlay = TRUE, names = NULL, ...){
  if(is.null(pars)){pars = names(x$stanfit)[1:4]
  message("no pars specified, plotting first four parameters")}

  np = length(pars)

  pars.x = intersect(pars, names(x$stanfit))
  if (length(pars.x) < 1) stop('no parameter found in model x')
  draws <- vector("list", np)
  names(draws) <- pars
  draws.raw = extract(x$stanfit, pars = pars.x)
  for (p in names(draws.raw)){
    draws[p] <- draws.raw[p]
  }
  if(!is.null(y)) {
    pars.y = intersect(pars, names(y$stanfit))
    if (length(pars.y) < 1) stop('no parameter found in model y')
    draws.y <- vector("list", np)
    names(draws.y)<- pars
    draws.raw.y = extract(y$stanfit, pars = pars.y)
    for (p in names(draws.raw.y)){
      draws.y[p] <- draws.raw.y[p]
    }
    }
  if(is.null(names))names = c('model 1','model 2')
  #cors<-round(cor(x$samples),2) #correlations

  # store old par
  old.par <- par(no.readonly=TRUE)

  # make layout for plot layout
  laymat<-diag(seq_len(np)) #histograms
  laymat[lower.tri(laymat)]<-(np + 1):(np + (np^2 - np) / 2) #correlations
  laymat[upper.tri(laymat)]<-(np + (np^2 - np) / 2 + 1):np^2 #heatmaps

  layout(laymat) #define layout using laymat

  par(mar=c(2,2,0.6,0.6), mgp = c(3,0.5,0)) #define marginals etc.

  # Draw histograms, tweak arguments of hist to make nicer figures
  for(i in seq_len(np)){
    if(is.null(y)){
      if (!is.null(draws[[i]])) { plot(density(draws[[i]]),main="", col = 'blue')
      } else {
        plot.default(type = 'n', ...)
        }
    } else {
      dens = switch(is.null(draws[[i]])+1,density(draws[[i]]),NULL)
      dens.y = switch(is.null(draws.y[[i]])+1,density(draws.y[[i]]),NULL)

      plot(dens ,main="", col = 'blue',
           xlim = range(dens$x,dens.y$x),
           ylim = range(dens$y,dens.y$y))
      lines(dens.y,main="", col = 'red')
      if(i == 1) legend('topleft', legend = names, text.col = c('blue','red'), bty = 'n')
      }
    title(main=pars[i], line = -0.1, xpd=TRUE)
  }

  #
  # Draw scatter
  for(i in seq_len(np-1))
    for(j in (i+1):np){
      if (scatter == TRUE){
        if(!is.null(draws[[i]]) & !is.null(draws[[j]])){

          plot.default(draws[[i]],draws[[j]], pch=16, cex=0.3, col='blue',
                       xlim = range(draws[[i]],draws.y[[i]]),
                       ylim = range(draws[[j]],draws.y[[j]]),
                       ...)
        } else {
          plot.default(range(draws[[i]]), range(draws[[j]]), type = 'n', ...)
          }
        if(!is.null(draws.y[[i]]) & !is.null(draws.y[[j]])){
          points(draws.y[[i]],draws.y[[j]], pch=16, cex=0.3, col='red', ...)}
      } else {
        plot.default(range(draws[[i]]), range(draws[[j]]), type = 'n', ...)
      }
    }

  # Plot heatmaps, here I use kde2d function for density estimation
  # image function for generating heatmaps
  #library(MASS)
  ## some pretty colors
  #library(RColorBrewer)
  k <- 5
  my.cols <- (RColorBrewer::brewer.pal(9, "Blues")[3:9])
  my.cols.y <- (RColorBrewer::brewer.pal(9, "Reds")[3:9])

  ## compute 2D kernel density, see MASS book, pp. 130-131


  for(i in 2:np)
    for(j in seq_len(i-1)){
      plot.default(range(draws[[i]],draws.y[[i]]), range(draws[[j]],draws.y[[j]]), type = 'n', ...)
      if(!is.null(draws[[i]]) & !is.null(draws[[j]])){
      z <- MASS::kde2d(draws[[i]],draws[[j]], n=25)
      contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE)}
      if(!is.null(y)){
        if(!is.null(draws.y[[i]]) & !is.null(draws.y[[j]])){
        z.y <- MASS::kde2d(draws.y[[i]],draws.y[[j]], n=20)
        contour(z.y, drawlabels=FALSE, nlevels=k, col=my.cols.y, add=TRUE)
      }}
      #if (medians) abline(h=median(x$samples[,j]), v=median(x$samples[,i]), lwd=2, lty = 2)
      #if (trend == TRUE) {
       # rows <- sample(nrow(x$samples), ceiling(0.05*nrow(x$samples)))
        #loess_fit <- loess(x$samples[rows,j] ~ x$samples[rows,i])
        #points(x$samples[rows,i], predict(loess_fit), col = "blue", pch=16, cex=0.5)
      #}

    }
  # restore old par
  par(old.par)
}
