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

##  PLAIN CLASSROOM EDITION
##  Generated from R/02_Visualization.R by scripts/make_class_src.py --
##  edit R/02_Visualization.R and regenerate; changes made here will be overwritten.
## ============================================================================

# R packages required
#install.packages("cluster")
library(cluster)

# Working directory
#   Point this at the folder that holds this session's data. Every file
#   name below is resolved relative to it, so the script and its data have
#   to travel together -- or at least stay in step.
#
#   This session reads: iris.data.csv  (optional)
setwd("YOUR/DIRECTORY")

# Built-in data sets used below
data(iris)

## ============================================================================
##  1. GETTING DATA INTO R
## ============================================================================
##  Where am I? Every path R sees is relative to this directory.
cat("Working directory:", getwd(), "\n")

##  The two workhorses are read.table() and read.csv(); read.csv() is just
##  read.table() with header = TRUE and sep = "," already set. Both take a path
##  relative to getwd(), which is why hard-coded absolute paths (the original
##  scripts began with setwd("C:\\Users\\...")) break on everyone else's machine.
csv <- "iris.data.csv"
if (file.exists(csv)) {
  m1 <- read.table(csv, header = TRUE, sep = ",")   # the explicit form
  m2 <- read.csv(csv)                               # the shorthand
  cat("read.table and read.csv agree:", identical(m1, m2), "\n")
  cat("columns in the course CSV:", paste(names(m2), collapse = ", "), "\n")
  ## Note: the course file puts the SPECIES in column 1 and names the
  ## measurements Sepal_L ... Petal_W. R's built-in iris puts Species LAST and
  ## uses Sepal.Length ... Petal.Width. Same data, different layout -- always
  ## look before you index by number.
}

##  Other formats you will meet: readRDS() for R's own binary format (that is
##  how leukemiaExpressionSubset.rds is stored), and read.csv(header = FALSE)
##  for a headerless numeric grid such as Limenitis_archippus.csv.

## ---- interrogating the object -----------------------------------------------
m <- iris

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
op <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
plot(X[, 1:2], pch = 19, col = "grey40", main = "No grouping")
plot(X[, 1:2], pch = 19, col = pal[grp],  main = "Coloured by species")
legend("topright", legend = levels(grp), fill = pal, bty = "n", cex = .8)
par(op)

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

par(mfrow = c(1, 1))
plot(X[, 1:2], type = "n", main = "Ellipsoid hulls and simulated Gaussians")
for (i in seq_along(hulls)) {
  points(sim_cloud(hulls[[i]]), col = pal[i], pch = ".")
  lines(predict(hulls[[i]]), col = pal[i], lwd = 2)
}
points(X[, 1:2], pch = 19, col = pal[grp], cex = .7)
legend("topright", legend = levels(grp), fill = pal, bty = "n", cex = .8)

## ---- 4. Marginal distributions ---------------------------------------------
## Bin width is a choice, not a property of the data: the same column looks
## unimodal at 8 bins and bimodal at 30. Say which you used.
op <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
hist(X[, 1], main = "Sepal.Length, default bins", xlab = "Sepal.Length",
     col = "grey85")
hist(X[, 1], nclass = 20, main = "Sepal.Length, 20 bins",
     xlab = "Sepal.Length", col = "grey85")
par(op)

## ---- 5. All pairs at once: the R view, panel by panel ----------------------
par(mfrow = c(1, 1))
pairs(X, pch = 19, col = pal[grp], cex = .7,
      main = "Iris: every pair of variables")

## ---- 5b. Three variables, rotatable ----------------------------------------
##  rgl opens a window you can drag with the mouse -- the best way to convey
##  that a 3-D scatter has no single correct viewpoint. It needs a live OpenGL
##  display, so it is skipped in batch runs and on machines without rgl.
if (interactive() && requireNamespace("rgl", quietly = TRUE)) {
  rgl::plot3d(X[, 1:3], col = pal[grp], type = "s", size = 1)
  rgl::rglwidget()          # embeds the widget when knitting; harmless otherwise
} else {
  message("Skipped: interactive 3-D scatter (drag to rotate).  Install with: install.packages(c(\"rgl\"))")
}

## ---- 6. Summaries by group, and exporting them -----------------------------
cat("\nColumn means (whole table):\n"); print(round(colMeans(X), 3))

out <- aggregate(X, by = list(Species = grp), FUN = mean)
cat("\nColumn means by species:\n"); print(out)

## Exporting: write.table() is the general form (you choose the separator),
## write.csv() the comma-specific shorthand. row.names = FALSE stops R from
## adding a leading column of "1", "2", "3" that no other program expects.
##
## Written to a temporary directory so sourcing this script never litters the
## repository; point `dest` at data/ or a results folder to keep the files.
dest_txt <- file.path(tempdir(), "means_iris.txt")
dest_csv <- file.path(tempdir(), "means_iris.csv")
utils::write.table(out, dest_txt, sep = "\t", row.names = FALSE)
utils::write.csv(out, dest_csv, row.names = FALSE)
cat("\nWrote", dest_txt, "\n      ", dest_csv, "\n")

cat("\n[02_Visualization] done.\n")
