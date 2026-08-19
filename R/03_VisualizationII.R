## ============================================================================
##  MS_LJMR :: 03_VisualizationII.R — Data Visualization II
##
##  Original authors: BIOL 943 course team
##                    (Jorge Soberón, Laura Jiménez & Marlon E. Cobos,
##                     University of Kansas)
##  Revised for the MS_LJMR edition by: Luis J. Madrigal-Roca
##
##  Revised from 03_VisualizationII.R ("Session 5: Visualization, Part 2").
##  The original set a working directory, read BiodivCountries.csv by hard-coded
##  COLUMN NUMBERS (3, 28:31, 34:36) and called x11(). Here the biodiversity
##  columns are selected by NAME — the same names 05_Clustering.R uses, which
##  the seeded fallback in 00_utils.R also emits — so the script survives any
##  reordering of the real file.
##
##  Where this sits in the course
##  -----------------------------
##  Session I showed the R view one pair at a time. Two variables fit on paper
##  and three fit in a projection; beyond that, plotting stops being an option
##  and ordination takes over. This session pushes direct display to its limit
##  so that the need for Chapters 2-3 is felt rather than asserted.
## ============================================================================
## ---------------------------------------------------------------------------
##  STANDALONE USE (no repository needed)
##  ---------------------------------------------------------------------------
##  As shipped, the source() line below borrows three helpers from the course
##  repository: get_data() (loads a data set), need() (loads packages) and
##  with_fig() (opens a plot device). To run this script entirely on its own,
##  delete that line and uncomment the block below. Nothing else changes.
##
##  The standalone with_fig() just draws each figure to the screen, one after
##  the other. To save them as files instead, replace its body with
##      png(paste0(name, ".png")); on.exit(dev.off()); force(expr)
##
##  Files to keep next to this script: BiodivCountries.csv
##
# need <- function(...) invisible(lapply(c(...), function(p) {
#   if (!requireNamespace(p, quietly = TRUE))
#     stop("This session needs: install.packages(\"", p, "\")", call. = FALSE)
#   suppressPackageStartupMessages(library(p, character.only = TRUE))
# }))
# with_fig <- function(name, expr, ...) invisible(force(expr))   # draw on screen
# get_data <- function(name) switch(name,
#   iris      = { utils::data("iris");   iris },
#   mtcars    = { utils::data("mtcars"); mtcars },
#   biodiv    = read.csv("BiodivCountries.csv", stringsAsFactors = TRUE),
#   stop("unknown data set: ", name))
## ---------------------------------------------------------------------------
if (!exists("get_data")) source(file.path("R", "00_utils.R"))
need("scatterplot3d", "car")

has3d <- requireNamespace("scatterplot3d", quietly = TRUE)
hascar <- requireNamespace("car", quietly = TRUE)

## ---- 1. Three dimensions, drawn honestly -----------------------------------
iris <- get_data("iris")
pal  <- c("#8c2d3a", "#1f3b73", "#2a7f7f")

if (has3d) {
  with_fig("03_iris_3d", {
    op <- par(mfrow = c(1, 2), mar = c(3, 3, 3, 1))
    scatterplot3d::scatterplot3d(iris[, 1:3], pch = 19, color = "grey40",
                                 main = "Three variables, no grouping")
    scatterplot3d::scatterplot3d(iris[, 1:3], pch = 19, color = pal[iris$Species],
                                 main = "Coloured by species")
    par(op)
  }, width = 11, height = 5.5)
} else {
  message("scatterplot3d not installed - skipping the 3-D iris panel.")
}

## A 3-D scatter is still a 2-D picture: the third axis is a projection chosen
## by `angle`, and the apparent structure changes with it. Rotating through a
## few angles is the cheapest demonstration that "seeing" >2 dimensions always
## means choosing a projection — which is exactly what PCA automates, except it
## chooses the projection by maximizing variance instead of by eye.
if (has3d) {
  with_fig("03_projection_angles", {
    op <- par(mfrow = c(1, 3), mar = c(3, 3, 3, 1))
    for (a in c(10, 40, 70))
      scatterplot3d::scatterplot3d(iris[, 1:3], pch = 19, color = pal[iris$Species],
                                   angle = a, main = paste("angle =", a))
    par(op)
  }, width = 13, height = 4.5)
}

## ---- 2. Aggregating a real table before plotting it ------------------------
## 186 countries is too many points to label. Aggregating to world regions
## turns an unreadable cloud into eight labelled points — the same move
## 05_Clustering.R makes before hierarchical clustering, and the same columns.
biodiv <- get_data("biodiv")
rich   <- c("AmphRich", "Rept_rich", "BirdRich", "MamsRich")
dens   <- c("DensAmphRich", "DensRept_rich", "DensBirdRich", "DensMamsRich")

## The real BiodivCountries.csv is incomplete, so aggregate with na.rm = TRUE
## rather than dropping whole countries for one missing value.
agg <- aggregate(biodiv[, c(rich, dens)],
                 by = list(Region = biodiv$RegionCode),
                 FUN = function(x) mean(x, na.rm = TRUE))
M   <- as.matrix(agg[, -1]); rownames(M) <- agg$Region
cat("Regional means:\n"); print(round(M[, rich], 1))

if (has3d) {
  with_fig("03_biodiv_3d", {
    s3d <- scatterplot3d::scatterplot3d(
      M[, c("AmphRich", "Rept_rich", "BirdRich")],
      pch = 19, color = "#1f3b73",
      type = "h",              # drop lines to the floor: restores depth cues
      xlab = "Amphibians", ylab = "Reptiles", zlab = "Birds",
      main = "Mean richness by world region", box = FALSE, angle = 25)
    ## xyz.convert() maps the 3-D coordinates onto the 2-D page, which is what
    ## makes it possible to add ordinary text() labels to a 3-D plot.
    co <- s3d$xyz.convert(M[, c("AmphRich", "Rept_rich", "BirdRich")])
    text(co$x, co$y, labels = rownames(M), cex = .8, pos = 3)
  })
}

## ---- 3. Scatterplot matrices with structure added --------------------------
## car::scatterplotMatrix() is pairs() plus a linear fit, a loess curve and a
## marginal density per panel. The loess is the useful part: it shows whether
## the linear summary that PCA and regression both assume is defensible.
cars <- get_data("mtcars")
if (hascar) {
  with_fig("03_spm_mtcars", {
    car::scatterplotMatrix(~ mpg + disp + drat + wt, data = cars,
                           regLine = list(col = "#8c2d3a"),
                           smooth  = list(col.smooth = "#1f3b73",
                                          col.spread = "#1f3b73"),
                           col = "grey30", pch = 19,
                           main = "mtcars: linear fit vs. loess")
  }, width = 8, height = 8)

  ## Conditioning on a factor splits every panel by group. When groups have
  ## different slopes, a single global ordination will average them away.
  with_fig("03_spm_by_cyl", {
    car::scatterplotMatrix(~ mpg + disp + drat + wt | factor(cyl), data = cars,
                           regLine = list(col = "grey20"), smooth = FALSE,
                           col = pal, pch = 19,
                           main = "mtcars, conditioned on cylinder count")
  }, width = 8, height = 8)
} else {
  message("car not installed - skipping the scatterplot-matrix panels.")
}

## ---- 4. Putting it together: labelled, coloured, conditioned 3-D -----------
if (has3d) {
  cyl_col <- pal[factor(cars$cyl)]
  with_fig("03_mtcars_3d", {
    s3d <- scatterplot3d::scatterplot3d(
      cars$disp, cars$wt, cars$mpg,
      color = cyl_col, pch = 19, type = "h", lty.hplot = 2, scale.y = .75,
      xlab = "Displacement (cu. in.)", ylab = "Weight (lb/1000)",
      zlab = "Miles per gallon", box = FALSE, angle = 20,
      main = "mtcars: three variables, conditioned on a fourth")
    co <- s3d$xyz.convert(cars$disp, cars$wt, cars$mpg)
    text(co$x, co$y, labels = rownames(cars), pos = 3, cex = .5)
    legend("topright", inset = .05, bty = "n", cex = .8,
           title = "Cylinders", legend = levels(factor(cars$cyl)), fill = pal)
  }, width = 9, height = 7)
}

## Four variables took every device available — colour, height, position and
## text. A fifth would not fit. That ceiling is the argument for ordination.
cat("\n[03_VisualizationII] done.\n")
