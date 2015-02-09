cols <- c (N = "#008000", `A°II` = "#0000FF", `A°III+` = "#FF0000", all = "black")
dimcols <- c (N = "#00800060", `A°II` = "#0000A060", `A°III+` = "#FF000060")

library (ggplot2)
library (hyperSpec)

load ("astrocytomas.RData")
astro


###### class membership values

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


memberships <- c (list (all = apply (astro$label, 1, max)),
                  apply (astro$label, 2, function (x) x [x > 0]))
before <- sapply (memberships, length)  # curiosity: how much reduction?
memberships <- lapply (memberships, fraction)
after <- sapply  (memberships, nrow)


memberships <- rbind.w.name (memberships)

summary (memberships)
dim (memberships)

png ("membership.png", width = 400, height = 250)
ggplot (data = memberships, aes(x = x, y = y, colour = class)) +
  geom_step () +
  scale_colour_manual ("class", value = cols) +
  scale_x_continuous (name = "Fraction of Spectra", expand = c (0, 0), limits = c (0, 1.005)) +
  scale_y_continuous (name = "Class Membership", expand = c (0, 0), limits = c (0, 1.005)) 
dev.off ()


##### spectra plot

library (Hmisc)

wtd.percentilespc <- function (weights, spc, probs = c (.16, .5, .84)){
  spc <- spc [weights > 0]              # reduce calculation time by removing unrelated spectra
  weights <- weights [weights > 0]

  spc <- apply (spc, 2, wtd.quantile, weights = weights, probs = probs, normwt = FALSE)
  spc$percentile <- probs

  spc
}

spc <- apply (astro$label, 2, wtd.percentilespc, spc = astro [, "spc"])
spc <- rbind.w.name (spc)

spc <- sweep (spc, 2, graunorm, `+`)

spc

df <- rbind (cbind (as.long.df (spc [,, min ~ 1800]), wlrange = "low"),
             cbind (as.long.df (spc [,, 2800 ~ max]), wlrange = "high"))

df <- cast (df, .wavelength + class + wlrange ~ percentile, value = "spc")


p <- ggplot (data = df) +
  geom_ribbon (aes (x = .wavelength, ymin = `0.16`, ymax = `0.84`, fill = class), col = "black", size = 0.15)  +
  scale_fill_manual ("class", value = dimcols) + 
  geom_line (aes (x = .wavelength, y = `0.5`), size = 0.25)

p <- p + facet_grid (class ~ wlrange, scales = "free_x", space = "free",
                     labeller = function (variable, value){
                       if (variable == "wlrange") "" else value
                     }) 

p <- p + ylab (labels (spc, "spc"))  +
  scale_x_continuous (name = labels (spc, ".wavelength"), breaks = seq (800, 3000, 200), expand = c (0, 50))  +
    opts (legend.position = "none")

png ("spc.png", width = 750, height = 300)
p
dev.off ()

####### 2d histogram of LDA projection

library (MASS)
lda <- lda (label.factor ~ spc, subset (astro$., label.factor != "soft"))
desc <- predict (lda, astro$.)$x


library (hexbin)

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

h <- hist2d (astro$label, desc, 75)

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

h$col <- I (sweep (h$counts, 2, apply (h$count, 2, max), `/`))
h$col <- colmix.rgb (h$col, purecol = cols [1:3])

p <- ggplot (data = h, aes (x = x, y = y, fill = col, colour = col, group = 1)) + geom_hex (stat = StatIdentity) +
   coord_equal () + scale_fill_identity() + scale_colour_identity() +
  opts (panel.background = theme_rect(fill = NA, colour = NA),
        panel.grid.major = theme_line(colour = NA), 
        panel.grid.minor = theme_line(colour = NA),
        plot.margin = unit(c(0.5, 0, 0 ,0), "lines")
        )  +
  scale_x_continuous (limits = c (-3.2,3.5)) +
  scale_y_continuous (limits = c (-3, 2)) +
  labs (x = "LDA 1", y = "LDA 2")

##### contours

library (sp)

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

median <- data.frame ()
for (class in seq_along (colnames (astro$label))){
  contours <- peel (desc, weights = astro$label [, class] > 1-1e-5, probs = c (0, .5))
  tmp <- apply (desc [contours [[1]],], 2, median)
  median <- rbind (median, data.frame (x = tmp [1], y = tmp [2], class = class))

  tmp <- data.frame (desc [contours [[2]],], class = class)
  colnames (tmp) <- c ("x", "y")
  p <- p + geom_polygon (data = tmp, aes (x = x, y = y), fill = NA, col = cols [class])  
}
p <- p + geom_point (data = median, aes (x = x, y = y), fill = "white", col = "black", shape = 21, size = 3)


####### legend

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


l <- legend (cols [-4], h$counts)

library (grid)
plot.with.legend.right <- function (graph, legend, legend.width = 8, legend.unit = "lines"){
  plot.new()
  pushViewport (viewport (layout = grid.layout (1, 2, 
                            widths = unit (c (1, legend.width), c("null",legend.unit))
                          )))
  print (graph,  viewport (layout.pos.col = 1), newpage = FALSE)
  print (legend, viewport (layout.pos.col = 2), newpage = FALSE)
  popViewport ()
}

png ("hist2d.png", width = 750, height = 450)
plot.with.legend.right (p, l)
dev.off ()


##### histogram hard spectra only
hard <- astro$label.factor != "soft"                        
h.hard <- hist2d(astro$label[hard,], desc[hard,], 75, rng = range (desc))

h.hard$col <- I (sweep (h.hard$counts, 2, apply (h$count, 2, max), `/`))
h.hard$col <- colmix.rgb (h.hard$col, purecol = cols [1:3])

p <- ggplot (data = h.hard, aes (x = x, y = y, fill = col, colour = col, group = 1)) + geom_hex (stat = StatIdentity) +
   coord_equal () + scale_fill_identity() + scale_colour_identity() +
  opts (panel.background = theme_rect(fill = NA, colour = NA),
        panel.grid.major = theme_line(colour = NA), 
        panel.grid.minor = theme_line(colour = NA),
        plot.margin = unit(c(0.5, 0, 0 ,0), "lines")
        )  +
  scale_x_continuous (limits = c (-3.2,3.5)) +
  scale_y_continuous (limits = c (-3, 2)) +
  labs (x = "LDA 1", y = "LDA 2")

median <- data.frame ()
for (class in seq_along (colnames (astro$label))){
  contours <- peel (desc, weights = astro$label [, class] > 1-1e-5, probs = c (0, .5))
  tmp <- apply (desc [contours [[1]],], 2, median)
  median <- rbind (median, data.frame (x = tmp [1], y = tmp [2], class = class))

  tmp <- data.frame (desc [contours [[2]],], class = class)
  colnames (tmp) <- c ("x", "y")
  p <- p + geom_polygon (data = tmp, aes (x = x, y = y), fill = NA, col = cols [class])  
}
p <- p + geom_point (data = median, aes (x = x, y = y), fill = "white", col = "black", shape = 21, size = 3)

png ("hist2dhard.png", width = 750, height = 450)
plot.with.legend.right (p, l)
dev.off ()

p
