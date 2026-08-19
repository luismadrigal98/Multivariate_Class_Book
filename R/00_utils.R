## ============================================================================
##  MS_LJMR :: 00_utils.R
##  Shared utilities for the Multivariate Data Analysis course code base
##  (Soberon / Jimenez / Cobos originals, revised & integrated with the ML
##   course by Luis J. Madrigal-Roca).
##
##  Design goals
##  ------------
##  * No hard-coded setwd() and no x11(): scripts run on Windows, macOS, Linux
##    and in batch (Rscript) sessions.
##  * Every data set is obtained through get_data(): it loads the real course
##    file when it is present on disk, and otherwise falls back to a
##    reproducible, seeded simulation with the SAME column structure, so the
##    scripts always run and figures always render. Drop the real CSV/RDS into
##    MS_LJMR/data/ (or set options(msljmr.data_dir=...)) and the loaders will
##    pick it up automatically -- no code change needed.
##  * open_dev() replaces x11(): interactive on screen, PNG in batch mode.
##
##  Author: Luis J. Madrigal-Roca   |   License: GPL-2 (matches course code)
## ============================================================================

## ---- paths -----------------------------------------------------------------
.msljmr_root <- function() {
  ## directory that holds this file's parent (…/MS_LJMR)
  getOption("msljmr.root", default = normalizePath(file.path(
    dirname(dirname(sys.frame(1)$ofile %||% ".")), "."), mustWork = FALSE))
}

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a[1])) b else a

data_dir <- function() {
  getOption("msljmr.data_dir",
            default = file.path(getOption("msljmr.root", "."), "data"))
}

fig_dir <- function() {
  d <- getOption("msljmr.fig_dir",
                 default = file.path(getOption("msljmr.root", "."), "figs"))
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

## ---- graphics device -------------------------------------------------------
##  open_dev("name", w, h): on an interactive session opens a screen device;
##  in batch (Rscript) writes figs/name.png. Call dev.off() when done, or use
##  with_fig().
open_dev <- function(name = "plot", width = 7, height = 6, res = 150) {
  if (interactive()) {
    grDevices::dev.new(width = width, height = height, noRStudioGD = TRUE)
  } else {
    grDevices::png(file.path(fig_dir(), paste0(name, ".png")),
                   width = width, height = height, units = "in", res = res)
  }
  invisible(NULL)
}

## convenience wrapper: with_fig("scree", { plot(...) })
with_fig <- function(name, expr, width = 7, height = 6, res = 150) {
  open_dev(name, width, height, res)
  on.exit(if (!interactive()) grDevices::dev.off())
  force(expr)
  invisible(NULL)
}

## ---- package helper --------------------------------------------------------
##  need("vegan","MASS"): load packages, telling the user how to install any
##  that are missing instead of crashing mid-script.
need <- function(...) {
  pkgs <- c(...)
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) {
    message("Missing packages: ", paste(miss, collapse = ", "),
            "\n  install.packages(c(",
            paste(sprintf('\"%s\"', miss), collapse = ", "), "))")
  }
  invisible(lapply(setdiff(pkgs, miss), function(p)
    suppressPackageStartupMessages(library(p, character.only = TRUE))))
}

## ---- evaluation helper -----------------------------------------------------
##  auc(p, y): area under the ROC curve from predicted scores `p` and binary
##  outcomes `y`, computed from the rank-sum identity AUC = (W - n1(n1+1)/2) /
##  (n1 n0) so no extra package is needed. Shared by the two ENM sessions.
auc <- function(p, y) {
  r  <- rank(p)
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

## ============================================================================
##  get_data(): the real-file-or-fallback loader
## ============================================================================
##  get_data("iris")                  -> Anderson's iris (built-in)
##  get_data("butterflies")           -> Quintana Roo wing-pattern counts
##  get_data("biodiv")                -> biodiversity-of-countries table
##  get_data("taxon")                 -> Crawley 4-taxa morphometrics
##  get_data("crawley")               -> speciesCrawley3 mixed-type table
##  get_data("europe")                -> 21 European cities (lon/lat)
##  get_data("countries_live")        -> ordinal liveability rankings
##  get_data("leukemia")              -> gene-expression subset (genes x samples)
##  get_data("limenitis")             -> grayscale wing image matrix
##  get_data("doubs")                 -> list(fish, env, xy) river gradient
## ----------------------------------------------------------------------------
get_data <- function(name, seed = 1998, quiet = FALSE) {
  name <- tolower(name)
  fmap <- list(
    butterflies    = "ButterfliesQRoo2.csv",
    biodiv         = "BiodivCountries.csv",
    taxon          = "taxon.csv",
    crawley        = "speciesCrawley3.csv",
    europe         = "CitiesEurope.csv",
    countries_live = "CountriesToLive.csv",
    leukemia       = "leukemiaExpressionSubset.rds",
    limenitis      = "Limenitis_archippus.csv",
    neotoma        = "NeotomaMorphoEnvir.csv",
    hanta          = "hanta_virtual.csv"
  )

  ## Datasets whose real file is a headerless numeric grid: no header row and no
  ## row names, every line pure data. read.csv()'s default header = TRUE would
  ## silently promote the first row of pixels to column names and hand back a
  ## data frame, which then fails in anything requiring a true matrix (norm(),
  ## svd() on a subset, image()). Read these with header = FALSE and coerce.
  matrix_files <- c("limenitis")
  ## 1) built-ins
  if (name == "iris")   { utils::data("iris",   envir = environment()); return(iris) }
  if (name == "mtcars") { utils::data("mtcars", envir = environment()); return(mtcars) }

  ## 2) real course file on disk?
  f <- fmap[[name]]
  if (!is.null(f)) {
    path <- file.path(data_dir(), f)
    if (file.exists(path)) {
      if (!quiet) message("get_data('", name, "'): using real file ", path)
      if (grepl("\\.rds$", f)) return(readRDS(path))
      if (name %in% matrix_files)
        return(as.matrix(utils::read.csv(path, header = FALSE)))
      return(utils::read.csv(path, stringsAsFactors = TRUE))
    }
  }

  ## 2a) GENUINE in-package datasets the class actually used: prefer the real
  ##     data when the package is installed on the user's machine.
  if (name == "doubs" && requireNamespace("ade4", quietly = TRUE)) {
    if (!quiet) message("get_data('doubs'): using the REAL ade4::doubs dataset")
    e <- new.env(); utils::data("doubs", package = "ade4", envir = e)
    return(list(fish = e$doubs$fish, env = e$doubs$env, xy = e$doubs$xy))
  }

  ## 2b) frozen DEMO file (data/demo/) -- shared byte-for-byte with the Python
  ##     verifier, so deterministic results (eigenvalues, inertias, loadings)
  ##     reproduce exactly. Generated by verify/freeze_and_reference.py.
  demo_dir  <- file.path(data_dir(), "demo")
  demo_read <- function(fn) utils::read.csv(file.path(demo_dir, fn), stringsAsFactors = TRUE)
  if (name == "doubs" && file.exists(file.path(demo_dir, "doubs_fish.csv"))) {
    if (!quiet) message("get_data('doubs'): using frozen demo data")
    return(list(fish = demo_read("doubs_fish.csv"), env = demo_read("doubs_env.csv")))
  }
  if (name == "limenitis" && file.exists(file.path(demo_dir, "limenitis.csv"))) {
    if (!quiet) message("get_data('limenitis'): using frozen demo data")
    return(as.matrix(utils::read.csv(file.path(demo_dir, "limenitis.csv"), header = FALSE)))
  }
  if (name == "leukemia" && file.exists(file.path(demo_dir, "leukemia.csv"))) {
    if (!quiet) message("get_data('leukemia'): using frozen demo data")
    m <- utils::read.csv(file.path(demo_dir, "leukemia.csv"), row.names = 1, check.names = FALSE)
    return(as.matrix(m))
  }
  demo1 <- c(taxon = "taxon.csv", crawley = "crawley.csv", hanta = "hanta.csv",
             pam = "pam.csv")
  if (!is.null(demo1[[name]]) && file.exists(file.path(demo_dir, demo1[[name]]))) {
    if (!quiet) message("get_data('", name, "'): using frozen demo data")
    return(demo_read(demo1[[name]]))
  }

  ## 3) reproducible fallback (only if no real or demo file is present)
  if (!quiet) message("get_data('", name,
                      "'): real file not found -- using seeded simulation.")
  set.seed(seed)
  sim <- switch(name,
    butterflies    = .sim_butterflies(),
    biodiv         = .sim_biodiv(),
    taxon          = .sim_taxon(),
    crawley        = .sim_crawley(),
    europe         = .sim_europe(),
    countries_live = .sim_countries_live(),
    leukemia       = .sim_leukemia(),
    limenitis      = .sim_limenitis(),
    doubs          = .sim_doubs(),
    neotoma        = .sim_neotoma(),
    hanta          = .sim_hanta(),
    pam            = .sim_pam(),
    stop("Unknown dataset: ", name))
  attr(sim, "simulated") <- TRUE
  sim
}

## ---- fallback generators (structure-faithful, seeded) ----------------------
##
##  A fallback is only useful if a script written against it also runs against
##  the real file, so each generator should emit the REAL file's column names.
##  Verified aligned: taxon, crawley, biodiv, hanta, limenitis, leukemia.
##  Known still divergent (no current script uses them, so they are documented
##  rather than changed blind -- align before writing a session that needs one):
##    butterflies    real: Genus, Sp, Name, Pattern, HD..PF
##                   sim : Superfamily, Genus, Species, Pattern, HD..PF
##    europe         real: City, Long, Lat          sim: City, lon, lat
##    countries_live real: Country, living, climate, food, security,
##                         hospitality, infrastructure
##                   sim : Country, Cost, Health, Safety, Climate, Freedom
##    neotoma        real: ID, sp, long, lat, No. cat, Code, d1..d8, ...
##                   sim : sp, lon, lat, m1..m80, bio1..bio19

## Quintana Roo butterflies: species x (meta + 7 successional sites)
.sim_butterflies <- function() {
  patterns <- c("Sand","BlackOrange","ShrubOrange","DarkCarpet","OpenContrast",
                "BlackRed","SilentBand","Bark","CanopyOrange","Adelpha",
                "CloseReflect","Tiger")
  sites <- c("HD","SD","GA","YA","MA","OA","PF")           # succession gradient
  ns <- 60
  species <- paste0("sp", sprintf("%02d", seq_len(ns)))
  pat <- sample(patterns, ns, replace = TRUE)
  ## a latent successional optimum per species drives a Poisson gradient
  opt <- runif(ns, 1, 7)
  M <- sapply(seq_along(sites), function(k)
    rpois(ns, lambda = 6 * exp(-0.5 * (k - opt)^2)))
  colnames(M) <- sites
  data.frame(Superfamily = "Nymphalidae", Genus = "Gen",
             Species = species, Pattern = pat, M, stringsAsFactors = TRUE)
}

## Biodiversity of countries: RegionCode + richness / area-corrected richness
## (density) / governance / capacity.
##
## Column names and RegionCode levels mirror the real BiodivCountries.csv so
## that a script written against the simulation runs unchanged on the real file.
## The real table has no per-taxon "diversity index"; what it carries alongside
## raw richness is Dens*Rich, richness per unit area. Earlier versions of this
## generator invented tidier *Div names that exist in no real file, which is
## exactly the mismatch that used to break 04_Clustering.R.
.sim_biodiv <- function() {
  regions <- c("APC","CAR","EECA","LAM","MENA","NAM","SSA","WEU")
  n <- 186
  reg <- sample(regions, n, replace = TRUE)
  base <- match(reg, regions)
  rich <- function(mult) round(abs(rnorm(n, 50 * mult + 8 * base, 25)))
  dens <- function(mult) round(abs(rnorm(n, 5 * mult + 0.8 * base, 3)), 3)
  gov  <- function() round(rnorm(n, 0, 1), 2)
  out <- data.frame(
    Country    = paste0("C", sprintf("%03d", seq_len(n))),
    RegionCode = factor(reg),
    AmphRich = rich(1.0), Rept_rich = rich(1.2),
    BirdRich = rich(3.0), MamsRich  = rich(1.5),
    DensAmphRich = dens(1.0), DensRept_rich = dens(1.2),
    DensBirdRich = dens(3.0), DensMamsRich  = dens(1.5),
    VoiceAccount = gov(), Stability = gov(), GovEffect = gov(), RuleLaw = gov(),
    Wealth = round(rlnorm(n, 9, 1)),
    Capacity = round(rnorm(n, 50 + 5 * base, 15)),
    stringsAsFactors = TRUE)
  ## the real file is not complete; imitate that so NA handling gets exercised
  for (j in c("AmphRich","Rept_rich","BirdRich","MamsRich",
              "DensAmphRich","DensRept_rich","DensBirdRich","DensMamsRich"))
    out[sample(n, 3), j] <- NA
  out
}

## Crawley "taxon": 4 taxa x 30, seven morphometric variables (Gaussian blobs)
.sim_taxon <- function() {
  g <- 4; per <- 30; p <- 7
  centers <- matrix(rnorm(g * p, 0, 3), g, p)
  X <- do.call(rbind, lapply(seq_len(g), function(k)
    matrix(rnorm(per * p), per, p) + matrix(centers[k, ], per, p, byrow = TRUE)))
  vn <- c("Petals","Internode","Sepal","Bract","Petiole","Leaf","Fruit")
  colnames(X) <- vn
  data.frame(Taxon = factor(rep(c("I","II","III","IV"), each = per)), X)
}

## Crawley "speciesCrawley3": counts + continuous env + soil factor
.sim_crawley <- function() {
  n <- 90
  Tmp    <- runif(n, -2, 12)
  Precip <- rgamma(n, 2, scale = 60)
  pH     <- runif(n, 3.5, 8.5)
  Soil   <- factor(sample(c("clay","loam","sand"), n, TRUE))
  mu_b    <- exp(1.6 + 0.04 * Tmp - 0.06 * (pH - 6)^2)   # positive mean
  Biomass <- rgamma(n, shape = 3, scale = mu_b / 3)       # right-skewed (non-normal)
  Species <- rpois(n, lambda = exp(1.2 + 0.06 * Tmp + 0.15 * (pH - 4)))
  data.frame(Species, Biomass, Tmp, Precip, pH, Soil)
}

## 21 European cities (real public coordinates) for PCoA / geographic MDS
.sim_europe <- function() {
  data.frame(
    City = c("Athens","Barcelona","Berlin","Brussels","Dublin","Geneva",
             "Helsinki","Lisbon","London","Madrid","Milan","Munich","Oslo",
             "Paris","Prague","Rome","Stockholm","Vienna","Warsaw","Zurich","Hamburg"),
    lon = c(23.73,2.17,13.40,4.35,-6.26,6.14,24.94,-9.14,-0.13,-3.70,9.19,11.58,
            10.75,2.35,14.42,12.50,18.07,16.37,21.01,8.54,9.99),
    lat = c(37.98,41.39,52.52,50.85,53.35,46.20,60.17,38.72,51.51,40.42,45.46,
            48.14,59.91,48.86,50.08,41.90,59.33,48.21,52.23,47.38,53.55))
}

## Ordinal liveability rankings (Q-mode ordinal example for NMDS)
.sim_countries_live <- function() {
  countries <- c("Norway","Switzerland","Canada","Germany","Japan","France",
                 "Spain","Portugal","Mexico","Brazil","India","Kenya","Vietnam")
  crit <- c("Cost","Health","Safety","Climate","Freedom")
  M <- sapply(crit, function(.) sample(seq_along(countries)))
  data.frame(Country = countries, M)
}

## Leukemia-like expression: genes x samples, 3 subtypes with block signatures
.sim_leukemia <- function() {
  types <- rep(c("ALL","AML","CLL"), each = 20)          # 60 samples
  g <- 300
  base <- matrix(rnorm(g * length(types), 6, 1), g, length(types))
  for (t in unique(types)) {
    idx <- which(types == t)
    sig <- sample(g, 40)
    base[sig, idx] <- base[sig, idx] + rnorm(1, 3, 0.3)
  }
  colnames(base) <- paste0(types, ".", ave(seq_along(types), types, FUN = seq_along))
  rownames(base) <- paste0("g", seq_len(g))
  base
}

## Synthetic grayscale "wing" image matrix (rank structure for SVD/NMF)
.sim_limenitis <- function(nr = 69, nc = 100) {
  x <- seq(-3, 3, length.out = nc); y <- seq(-2, 2, length.out = nr)
  G <- outer(y, x, function(yy, xx)
    exp(-(xx^2)/3 - (yy^2)/1.2) * (1 + 0.6 * cos(3 * xx)) )
  G <- G + 0.15 * outer(y, x, function(yy, xx) sin(2 * xx) * cos(2 * yy))
  G <- (G - min(G)) / (max(G) - min(G))
  round(G, 4)
}

## Doubs-like river gradient: fish counts, env, spatial coords
.sim_doubs <- function() {
  n <- 30; s <- 27
  grad <- seq(0, 1, length.out = n)             # up- to down-stream
  opt  <- seq(0, 1, length.out = s)
  fish <- sapply(seq_len(s), function(k)
    rpois(n, lambda = 8 * exp(-((grad - opt[k])^2) / 0.02)))
  colnames(fish) <- paste0("Sp", seq_len(s))
  env <- data.frame(
    alt = 900 - 850 * grad + rnorm(n, 0, 10),
    slo = exp(4 - 4 * grad + rnorm(n, 0, .3)),
    flo = 10 + 400 * grad + rnorm(n, 0, 20),
    pH  = 8 - grad + rnorm(n, 0, .1),
    oxy = 12 - 4 * grad + rnorm(n, 0, .5),
    nit = 1 + 3 * grad + rnorm(n, 0, .2))
  xy <- data.frame(x = cumsum(runif(n, .5, 1.5)), y = cumsum(rnorm(n, 0, .4)))
  list(fish = as.data.frame(fish), env = env, xy = xy)
}

## Neotoma-like morphology + bioclim (for PCA parts/procrustes)
.sim_neotoma <- function() {
  n <- 615; sp <- factor(sample(paste0("N", 1:15), n, TRUE))
  base <- as.integer(sp)
  morph <- sapply(1:80, function(j) rnorm(n, 10 + 0.2 * base, 2))
  bio <- sapply(1:19, function(j) rnorm(n, 100 + 5 * base, 30))
  colnames(morph) <- paste0("m", 1:80); colnames(bio) <- paste0("bio", 1:19)
  data.frame(sp = sp, lon = runif(n, -115, -70), lat = runif(n, 20, 45),
             morph, bio)
}

## Hantavirus-like presence/absence for ENM (binomial on two bioclim predictors)
.sim_hanta <- function() {
  n <- 400
  bio_1  <- runif(n, -5, 30); bio_12 <- runif(n, 100, 2500)
  eta <- -2 + 0.15 * bio_1 - 0.0015 * bio_12 - 0.004 * (bio_1 - 15)^2
  Sp  <- rbinom(n, 1, 1 / (1 + exp(-eta)))
  data.frame(Sp = Sp, bio_1 = bio_1, bio_12 = bio_12)
}

## Presence-absence matrix: sites x species along a latitudinal gradient
.sim_pam <- function() {
  nsite <- 80; nsp <- 120
  lat <- seq(0, 1, length.out = nsite)
  opt <- runif(nsp); br <- runif(nsp, 0.05, 0.25)
  prob <- exp(-outer(lat, opt, `-`)^2 / (2 * rep(br, each = nsite)^2))
  pam <- matrix(as.integer(runif(nsite * nsp) < prob), nsite, nsp)
  colnames(pam) <- paste0("sp", seq_len(nsp))
  data.frame(lat = lat, pam)
}

message("MS_LJMR utils loaded. get_data(<name>), open_dev(), with_fig(), need().")
