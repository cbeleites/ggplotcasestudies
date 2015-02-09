library (ggplot2)
library (hyperSpec)
library (grid)
library (Hmisc)
library (MASS)
library (hexbin)
library (sp)

source ("ggplot-gliomas-functions.R")
load ("astrocytomas.RData")
load ("astro.RData")
astro

cols <- c (N = "#008000", `A°II` = "#0000FF", `A°III+` = "#FF0000", all = "black")
dimcols <- c (N = "#00800060", `A°II` = "#0000A060", `A°III+` = "#FF000060")

###### class membership values

memberships <- c (list (all = apply (astro$label, 1, max)),
                  apply (astro$label, 2, function (x) x [x > 0]))
before <- sapply (memberships, length)  # curiosity: how much reduction?
memberships <- lapply (memberships, fraction)
after <- sapply  (memberships, nrow)


memberships <- rbind.w.name (memberships)

summary (memberships)
dim (memberships)

png ("membership.png", width = 600, height = 350, res = 100)
ggplot (data = memberships, aes(x = x, y = y, colour = class)) +
  geom_step () +
  scale_colour_manual ("class", value = cols) +
  scale_x_continuous (name = "Fraction of Spectra", expand = c (0, 0), limits = c (0, 1.005)) +
  scale_y_continuous (name = "Class Membership", expand = c (0, 0), limits = c (0, 1.005)) 
dev.off ()


##### spectra plot

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

png ("spc.png", width = 750, height = 450, res = 100)
p
dev.off ()

####### 2d histogram of LDA projection

load ("astrocytomas.RData")
lda <- lda (label.factor ~ spc, subset (astro$., label.factor != "soft"))
desc <- predict (lda, astro$.)$x

h <- hist2d (astro$label, desc, 75)
h$col <- I (sweep (h$counts, 2, apply (h$count, 2, max), `/`))
h$col <- colmix.rgb (h$col, purecol = cols [1:3])

p <- ggplot (data = h, aes (x = x, y = y, fill = col, colour = col, group = 1)) + geom_hex (stat = StatIdentity) +
   coord_equal () + scale_fill_identity() + scale_colour_identity() +
  opts (panel.background = theme_rect(fill = NA, colour = "gray75"),
        panel.grid.major = theme_line(colour = NA), 
        panel.grid.minor = theme_line(colour = NA),
        plot.margin = unit(c(0.5, 0, 0 ,0), "lines")
        )  +
  scale_x_continuous (limits = c (-3.2,3.5)) +
  scale_y_continuous (limits = c (-3, 2)) +
  labs (x = "LDA 1", y = "LDA 2")

##### contours

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

l <- legend (cols [-4], h$counts)

png ("hist2d.png", width = 750, height = 450, res = 100)
plot.with.legend.right (p, l)
dev.off ()


##### histogram hard spectra only
hard <- astro$label.factor != "soft"                        
h.hard <- hist2d(astro$label[hard,], desc[hard,], 75, rng = range (desc))

h.hard$col <- I (sweep (h.hard$counts, 2, apply (h$count, 2, max), `/`))
h.hard$col <- colmix.rgb (h.hard$col, purecol = cols [1:3])

p <- ggplot (data = h.hard, aes (x = x, y = y, fill = col, colour = col, group = 1)) + geom_hex (stat = StatIdentity) +
   coord_equal () + scale_fill_identity() + scale_colour_identity() +
  opts (panel.background = theme_rect(fill = NA, colour = "gray75"),
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

png ("hist2dhard.png", width = 750, height = 450, res = 100)
plot.with.legend.right (p, l)
dev.off ()


