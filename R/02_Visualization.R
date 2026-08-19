## ============================================================================
##  MS_LJMR :: 02_Visualization.R — R fundamentals & Data Visualization I
##
##  Original authors: Laura Jiménez & Jorge Soberón
##                    (BIOL 943 Multivariate Methods, University of Kansas)
##  Revised for the MS_LJMR edition by: Luis J. Madrigal-Roca
##
##  Revised from 02_Visualization.R / 02_VisualizationB.R ("Swimming in the R
##  pool"). The original opened with setwd() to a personal Dropbox folder and
##  read iris.data.csv from it, then called x11() four times; both are gone.
##  Data comes through get_data(), figures through with_fig(), and the exported
##  table is written to a session temporary directory instead of the user's
##  working directory.
##
##  Why this session is in a multivariate course
##  --------------------------------------------
##  Everything later in the book is a projection of a data matrix. Before
##  projecting anything, you have to be able to READ the matrix: how many rows
##  and columns, which columns are numeric and which are factors, what the
##  marginal distributions look like, and how the variables covary in pairs.
##  The scatterplot matrix at the end of this script is already the "R view" of
##  Chapter 1 drawn one panel at a time.
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
##  Files to keep next to this script: none — the data sets used here ship with R itself
##
# need <- function(...) invisible(lapply(c(...), function(p) {
#   if (!requireNamespace(p, quietly = TRUE))
#     stop("This session needs: install.packages(\"", p, "\")", call. = FALSE)
#   suppressPackageStartupMessages(library(p, character.only = TRUE))
# }))
# with_fig <- function(name, expr, ...) invisible(force(expr))   # draw on screen
# get_data <- function(name) { utils::data("iris");   iris }
## ---------------------------------------------------------------------------
if (!exists("get_data")) source(file.path("R", "00_utils.R"))
need("cluster")

## ---- 1. Reading and interrogating a data matrix ----------------------------
m <- get_data("iris")

cat("mode:", mode(m), " class:", class(m), "\n")
cat("dim :", nrow(m), "rows x", ncol(m), "columns\n\n")
print(head(m))
str(m)
print(summary(m))
cat("\nGroups in the last column:\n"); print(unique(m$Species))

X   <- m[, 1:4]                                  # the numeric matrix
grp <- factor(m$Species)                         # the grouping factor
pal <- c("#8c2d3a", "#1f3b73", "#2a7f7f")        # one colour per species

## ---- 2. The scatterplot, plain and grouped ---------------------------------
## The plain version answers "is there structure?"; the grouped version answers
## "is the structure the one I already know about?" — the difference between
## unsupervised and supervised thinking, in two lines of code.
with_fig("02_scatter_grouped", {
  op <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  plot(X[, 1:2], pch = 19, col = "grey40", main = "No grouping")
  plot(X[, 1:2], pch = 19, col = pal[grp],  main = "Coloured by species")
  legend("topright", legend = levels(grp), fill = pal, bty = "n", cex = .8)
  par(op)
}, width = 10, height = 5)

## ---- 3. Ellipsoid hulls: the first multivariate summary --------------------
## A univariate summary is a mean and an SD. Its multivariate counterpart is a
## centroid and a covariance matrix, and the natural picture of that pair is an
## ellipse. cluster::ellipsoidhull() fits the minimum-volume ellipsoid to each
## species, giving centre ($loc) and shape ($cov).
hulls <- lapply(levels(grp), function(g)
  cluster::ellipsoidhull(as.matrix(X[grp == g, 1:2])))
names(hulls) <- levels(grp)

## The same covariance matrix also GENERATES data, via the affine property of
## the multivariate normal: if Sigma = A'A (Cholesky) and z ~ N(0, I), then
## A'z + mu ~ N(mu, Sigma). Overlaying such a cloud on the real points is a
## quick visual test of how well a Gaussian describes each species.
sim_cloud <- function(h, n = 4000) {
  A <- chol(h$cov)
  matrix(rnorm(2 * n), n, 2) %*% A +
    matrix(h$loc, n, 2, byrow = TRUE)
}
set.seed(1998)

with_fig("02_ellipsoid_hulls", {
  plot(X[, 1:2], type = "n", main = "Ellipsoid hulls and simulated Gaussians")
  for (i in seq_along(hulls)) {
    points(sim_cloud(hulls[[i]]), col = pal[i], pch = ".")
    lines(predict(hulls[[i]]), col = pal[i], lwd = 2)
  }
  points(X[, 1:2], pch = 19, col = pal[grp], cex = .7)
  legend("topright", legend = levels(grp), fill = pal, bty = "n", cex = .8)
})

## ---- 4. Marginal distributions ---------------------------------------------
## Bin width is a choice, not a property of the data: the same column looks
## unimodal at 8 bins and bimodal at 30. Say which you used.
with_fig("02_histograms", {
  op <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  hist(X[, 1], main = "Sepal.Length, default bins", xlab = "Sepal.Length",
       col = "grey85")
  hist(X[, 1], nclass = 20, main = "Sepal.Length, 20 bins",
       xlab = "Sepal.Length", col = "grey85")
  par(op)
}, width = 10, height = 5)

## ---- 5. All pairs at once: the R view, panel by panel ----------------------
with_fig("02_pairs", {
  pairs(X, pch = 19, col = pal[grp], cex = .7,
        main = "Iris: every pair of variables")
})

## ---- 6. Summaries by group, and exporting them -----------------------------
cat("\nColumn means (whole table):\n"); print(round(colMeans(X), 3))

out <- aggregate(X, by = list(Species = grp), FUN = mean)
cat("\nColumn means by species:\n"); print(out)

## Written to a temporary directory so sourcing this script never litters the
## repository; point `dest` at data/ or a results folder to keep the file.
dest <- file.path(tempdir(), "means_iris.csv")
utils::write.csv(out, dest, row.names = FALSE)
cat("\nWrote", dest, "\n")

cat("\n[02_Visualization] done.\n")
