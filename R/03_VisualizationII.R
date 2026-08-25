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
##  columns are selected by NAME — the same names 04_Clustering.R uses, which
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
##  As shipped, the source() line below borrows a few helpers from the course
##  repository: get_data() (loads a data set), need() (loads packages),
##  with_fig() (opens a plot device) and has_pkg()/skip_note() (let an optional
##  section be skipped rather than crash). To run this script entirely on its
##  own, delete that line and uncomment the block below. Nothing else changes.
##
##  The standalone with_fig() draws each figure to the current device -- inside
##  RStudio that is the Plot pane, so figures accumulate in the plot history.
##  To save them as files instead, replace its body with
##      png(paste0(name, ".png")); on.exit(dev.off()); force(expr)
##
##  Files to keep next to this script: BiodivCountries.csv; BiodiversityCountriesPCValues.csv
##
# need <- function(...) invisible(lapply(c(...), function(p) {
#   if (!requireNamespace(p, quietly = TRUE))
#     stop("This session needs: install.packages(\"", p, "\")", call. = FALSE)
#   suppressPackageStartupMessages(library(p, character.only = TRUE))
# }))
# with_fig <- function(name, expr, ...) {          # draw on the current device;
#   op <- par(no.readonly = TRUE); on.exit(par(op))  # in RStudio that is the Plot pane
#   invisible(force(expr))
# }
# has_pkg <- function(...) all(vapply(c(...), requireNamespace, logical(1), quietly = TRUE))
# skip_note <- function(what, pkgs) {
#   message("  [skipped] ", what, " -- needs ", paste(pkgs, collapse = ", "),
#           ":  install.packages(c(",
#           paste(sprintf('"%s"', pkgs), collapse = ", "), "))")
#   invisible(FALSE)
# }
# get_data <- function(name) switch(name,
#   biodiv    = read.csv("BiodivCountries.csv", stringsAsFactors = TRUE),
#   biodiv_pc = read.csv("BiodiversityCountriesPCValues.csv", stringsAsFactors = TRUE),
#   iris      = { utils::data("iris");   iris },
#   mtcars    = { utils::data("mtcars"); mtcars },
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
    par(mfrow = c(1, 2), mar = c(3, 3, 3, 1))
    scatterplot3d::scatterplot3d(iris[, 1:3], pch = 19, color = "grey40",
                                 main = "Three variables, no grouping")
    scatterplot3d::scatterplot3d(iris[, 1:3], pch = 19, color = pal[iris$Species],
                                 main = "Coloured by species")
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
    par(mfrow = c(1, 3), mar = c(3, 3, 3, 1))
    for (a in c(10, 40, 70))
      scatterplot3d::scatterplot3d(iris[, 1:3], pch = 19, color = pal[iris$Species],
                                   angle = a, main = paste("angle =", a))
  }, width = 13, height = 4.5)
}

## ---- 2. Aggregating a real table before plotting it ------------------------
## 186 countries is too many points to label. Aggregating to world regions
## turns an unreadable cloud into eight labelled points — the same move
## 04_Clustering.R makes before hierarchical clustering, and the same columns.
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

## ============================================================================
##  5. THE GRAMMAR OF GRAPHICS
## ============================================================================
##  Everything above is base graphics: you issue drawing commands in order.
##  ggplot2 works differently -- you declare a MAPPING from variables to visual
##  channels (aes) and then add layers (geoms) that render it. The payoff is
##  that conditioning on a factor stops being a loop and becomes one more term.
##
##  Data: the biodiversity table reduced to its principal components, which is
##  also a preview of Chapter 4 -- Biodiversity is PC1 of the area-corrected
##  richness variables (75.9% of their variance), Governance summarises 14
##  variables (86.7%), Wealth 3 (55.7%) and Capacity_no_GEF 19 (55.3%).
if (has_pkg("ggplot2")) {
  library(ggplot2)
  m3 <- get_data("biodiv_pc")
  cat("\nBiodiversity PC table:", nrow(m3), "countries x", ncol(m3), "columns\n")
  cat("  ", paste(names(m3), collapse = ", "), "\n")

  ## ggplot() returns an OBJECT. At the console it prints itself and you see the
  ## plot; inside a function or a loop you must print() it explicitly, or
  ## nothing is drawn. This catches everyone once.
  p1 <- ggplot(m3, aes(x = Biodiversity))

  with_fig("03_gg_histogram", {
    print(p1 + geom_histogram(colour = "black", fill = "#e8c33a", bins = 30) +
            theme_bw(base_size = 14) +
            labs(title = "One geom, one variable"))
  })

  ## the same mapping, a different geom: density instead of counts
  with_fig("03_gg_density", {
    print(p1 + geom_density(colour = "black", linewidth = 1.2, fill = "#2a7f7f") +
            theme_bw(base_size = 14) +
            labs(title = "geom_density: the smoothed version"))
  })

  ## ---- boxplots against a categorical variable -----------------------------
  ##  Hinges are the 1st and 3rd quartiles; whiskers reach the furthest point
  ##  within 1.5 x IQR of the hinge; anything past them is drawn as a point
  ##  (Tukey). The notch spans 1.58 x IQR/sqrt(n) and gives a rough 95% interval
  ##  for the median, so non-overlapping notches suggest different medians
  ##  (McGill et al. 1978).
  with_fig("03_gg_boxplots", {
    print(ggplot(m3, aes(x = RegionCode, y = Biodiversity)) +
            geom_boxplot(outlier.shape = 3, notch = TRUE) +
            theme_bw(base_size = 13) +
            labs(title = "Biodiversity by world region"))
  }, width = 9, height = 5)

  ## the same plot with the raw points jittered on top, coloured by region:
  ## a boxplot hides n, and n is often the whole story
  with_fig("03_gg_boxplot_jitter", {
    print(ggplot(m3, aes(x = GDPgroup, y = Biodiversity)) +
            geom_boxplot(outlier.shape = 3, notch = TRUE) +
            geom_jitter(aes(colour = RegionCode),
                        position = position_jitter(width = .2), size = 2) +
            theme_bw(base_size = 13) +
            labs(title = "Boxplot + the points it summarises"))
  }, width = 9, height = 5)

  ## ---- facets: the same plot, once per group -------------------------------
  ##  facet_wrap() is the grammar's answer to "does this relationship hold in
  ##  every region?" -- the question a single pooled regression cannot answer.
  with_fig("03_gg_facets", {
    print(ggplot(m3, aes(x = Wealth, y = Biodiversity)) +
            geom_point(aes(colour = GDPgroup)) +
            stat_smooth(method = lm, formula = y ~ x,
                        colour = "#8c2d3a", fill = "#e8c33a") +
            facet_wrap(~ RegionCode) + theme_bw(base_size = 11) +
            labs(title = "Wealth vs biodiversity, one panel per region"))
  }, width = 10, height = 7)

  ## facet_grid() crosses two factors instead of wrapping one
  with_fig("03_gg_facet_grid", {
    print(ggplot(m3, aes(x = Governance, y = Biodiversity)) +
            geom_point(aes(colour = RegionCode)) +
            stat_smooth(method = lm, formula = y ~ x, colour = "#8c2d3a") +
            facet_grid(Landlock ~ GDPgroup) + theme_bw(base_size = 11) +
            labs(title = "Landlocked (rows) x GDP group (columns)"))
  }, width = 10, height = 7)

  ## ---- a continuous variable mapped to colour ------------------------------
  with_fig("03_gg_gradient", {
    print(ggplot(m3, aes(x = Governance, y = Capacity_no_GEF)) +
            geom_point(aes(colour = Biodiversity), size = 3) +
            geom_point(shape = 1, size = 3, colour = "grey30") +
            scale_colour_gradientn(colours = hcl.colors(5, "Viridis")) +
            theme_bw(base_size = 13) +
            labs(title = "Three variables, no third axis"))
  })
} else {
  skip_note("the ggplot2 section", "ggplot2")
}

## Four variables took every device available — colour, height, position and
## text. A fifth would not fit. That ceiling is the argument for ordination.
cat("\n[03_VisualizationII] done.\n")
