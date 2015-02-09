### Helper functions for the ggplot2 case study
###
### (c) 2010 Claudia Beleites, cbeleites at units dot it
###
### This code is released under GPL 3 http://www.gnu.org/licenses/gpl-3.0.html
###

## function to tabulate values against fraction of objects
fraction <- function (max.membership){
  max.membership <- round (max.membership, digits = 3)
  y <- sort (unique (max.membership), decreasing = TRUE)
  x <- tabulate (factor (max.membership, levels = y))
  x <- cumsum (x)
  x <- x / tail (x, 1)
  x <- c (0, x)
  y <- c (y [1], y)

  data.frame (x = x, y = y)
}

# function to rbind objects in a list
# putting the names of the list elements as grouping factor into column 
rbind.w.name <- function (l){
  n <- names (l)
  for (i in n)
    l [[i]]$class <- i

  l <- do.call (rbind, l)
  l$class <- factor (l$class, levels = n)

  l
} 

# helper function for the weighted percentile spectra
wtd.percentilespc <- function (weights, spc, probs = c (.16, .5, .84)){
  spc <- spc [weights > 0]              # reduce calculation time by removing unrelated spectra
  weights <- weights [weights > 0]

  spc <- apply (spc, 2, wtd.quantile, weights = weights, probs = probs, normwt = FALSE)
  spc$percentile <- probs

  spc
}


# colour mixing
colmix.rgb <- function (x, purecol, against = 1, min = 0, max = 1, sub = TRUE){
  if (is.character (purecol))
    purecol <- t (col2rgb (purecol)) / 255

  if (sub)
    x <- against - x %*% (against - purecol)
  else
    x <- x %*% purecol
  
  x [x < min] <- min
  x [x > max] <- max

  cols <- rep (NA, nrow (x))
  cols [! is.na (x [,1])] <-   rgb (x [!is.na (x [, 1]),])

  cols
}

# 2d histograms for multiclass data
hist2d <- function (labels, desc, xbins, rng = range (desc)){
  h <- hexbin (desc, xbins = xbins,
               xbnds = rng, ybnds = rng,
               IDs = TRUE)
  
  counts <- hexTapply (h, 1 : nrow (labels),
                       function (i, labels){
                         colSums (labels [i,, drop = FALSE])
                       },
                       labels = labels)
  counts <- t (matrix (unlist (counts), nrow = length (counts [[1]])))
  colnames (counts) <- colnames (labels)

  data.frame (hcell2xy (h), counts = I(counts))
}

# legend for the multiclass histogram
legend <- function (purecolours, counts, dx = 0.33, classlabels = names (purecolours)) {
  if (! is.matrix (counts))
    counts <- matrix (counts, ncol = length (counts))
  
  maxcnt <- apply (counts, 2, max)
  df <- data.frame ()
  for (class in seq_along (maxcnt)){
    if (max (maxcnt) == 1)
      tmp <- c (0, seq_len (maxcnt [class] * 100) / 100)
    else
      tmp <- c (0, seq_len (maxcnt [class]))
    
    df <- rbind (df, data.frame (class = class,
                                 col = colmix.rgb (tmp / maxcnt [class], purecolours [class]),
                                 counts = tmp,
                                 dx = dx,
                                 dy = tmp[2] - tmp [1]))
  }


  l <- ggplot (df, aes (x = class), col = col) +
    geom_point (aes (x=class, y = 1), col = NA)  # trick to access continuous x values
  l <- l + geom_rect (aes (xmin = as.numeric (class) - dx,
                           xmax = as.numeric (class) + dx,
                           ymin = counts - dy,
                           ymax = counts + dy,
                           fill = col,
                           colour = col), dx = force (dx)
                      )

  l <- l + opts (plot.margin = unit(c(0.5, 0, 0 ,0), "lines"),
                 axis.text = p$options$legend.text,
                 axis.title = p$options$legend.text 
                 ) +
       scale_fill_identity () + scale_colour_identity () +
       scale_y_continuous (name = "counts", expand = c(0, max (df$dy)),
                           minor = pretty (c (0, max (maxcnt)), 25)) +
       scale_x_continuous (name = "class", expand = c (0, 0.5 * dx), minor = NA, major = NA,
                           breaks = seq_along (maxcnt), labels = classlabels)

  l
}

# peeling contours: 2d quantiles
peel <- function (x, y, weights = 1, probs = NA, threshold = 1 - 1e-3){
  if (missing (y) && ncol (x) == 2){
    y <- x [, 2]
    x <- x [, 1]
  }

  if  (length (x) != length (y))
    stop ("x and y need to have the same length.")

  weights <- rep (weights, length.out = length (x))
  
  ## start with all points
  pts.in <- seq_along (x)
  step <- 1
  hulls <- list ()

  ## too small weights can confuse the peeling as the hull polygon treats all points equally

  exclude <- weights < threshold
  if (any (exclude)) {
  ##   warning (sum (exclude), " points put into first hull due to too small weights")

    hulls [[1]] <- pts.in [exclude]
    pts.in <- pts.in [! exclude]
    step <- step + 1
  }
    
  ## peel off the hull polygons until nothing is left
  while (length (pts.in) > 1){
    hull <- chull (x [pts.in], y [pts.in])
    hulls [[step]] <- pts.in [hull]
    pts.in <- pts.in [-hull]
    step <- step + 1
  }

  # now count the number of point-equivalents in each hull
  n <- sapply (hulls, function (i) sum (weights [i]))

  ## and convert to percentiles
  n <- cumsum (n)
  qtl <- c(1, 1 - head (n, -1) / tail (n, 1))
  
  names (hulls) <- qtl
  
  if (! all (is.na (probs))){
    i <- round (approx (qtl[-1], seq_along (hulls[-1]), probs, rule = 2)$y) + 1
    hulls <- hulls [i]
  }
  
  hulls
}

# plotting with custom legend on the right side
plot.with.legend.right <- function (graph, legend, legend.width = 8, legend.unit = "lines"){
  plot.new()
  pushViewport (viewport (layout = grid.layout (1, 2, 
                            widths = unit (c (1, legend.width), c("null",legend.unit))
                          )))
  print (graph,  viewport (layout.pos.col = 1), newpage = FALSE)
  print (legend, viewport (layout.pos.col = 2), newpage = FALSE)
  popViewport ()
}
